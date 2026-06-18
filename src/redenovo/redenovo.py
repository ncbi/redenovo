import logging
import sys
import os
import pandas as pd
import numpy as np
import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
import time
from scipy.sparse.csgraph import connected_components
import multiprocessing as mp

from redenovo import __version__
from redenovo import arguments
from redenovo import resources
from redenovo import optimizer as model
from redenovo import exporter

__author__ = "ReDeNovo"
__copyright__ = "ReDeNovo"
__license__ = "GPL-3.0-only"
_logger = logging.getLogger(__name__)

_WORKER_STATE = {}

tool_dir = os.path.dirname(os.path.abspath(__file__))

def setup_logging(verbosity):
    """Setup basic logging."""

    logformat = "[%(asctime)s] %(levelname)s:%(name)s: %(message)s"

    if isinstance(verbosity, str):
        level = getattr(logging, verbosity.upper(), None)

        if not isinstance(level, int):
            raise ValueError(f"Invalid logging verbosity: {verbosity}")

    elif isinstance(verbosity, int):
        level = verbosity

    else:
        raise TypeError(f"Invalid logging verbosity type: {type(verbosity)}")

    logging.basicConfig(
        level=level,
        stream=sys.stdout,
        format=logformat,
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    for handler in root_logger.handlers:
        handler.setLevel(level)

    # Prevent TensorFlow messages from being printed twice.
    tf_logger = logging.getLogger("tensorflow")
    tf_logger.handlers.clear()
    tf_logger.propagate = True
    tf_logger.setLevel(logging.WARNING)
    

def main(args):
    print(f"ReDeNovo version: {__version__}")
    args = arguments.parse_args(args)
    setup_logging(args.verbosity)

    print_config(args)
    os.makedirs(args.out, exist_ok=True)
    stoppingCriteria = False
    
    start_time = time.time()
    current_numpri = args.numpri
    current_args_out = args.out
    novel_count = 0
                    
    novel_signatures = None
    if args.manual_cosmic:
        try:
            sigDB = pd.read_table(args.manual_cosmic_file, header=0)
            COSMIC = sigDB.iloc[:, 1:].T.apply(pd.to_numeric)
            COSMIC.index = sigDB.columns[1:]
            COSMIC.columns = sigDB.iloc[:, 0].astype(str)
            
            if len(args.exclude) > 0:
                missing_exclude = [sig for sig in args.exclude if sig not in COSMIC.index]
                if missing_exclude:
                    raise ValueError(f"Excluded signature(s) not found in manual COSMIC: {missing_exclude}")
                COSMIC = COSMIC.drop(args.exclude)
                
            print("Manual COSMIC loaded:", COSMIC.shape)

        except Exception as e:
            raise RuntimeError(f"Failed to load manual COSMIC file: {e}")

    else:
        COSMIC = get_COSMIC(args.genome, args.whole, args.exclude, args.cosmic_version)

    if args.add_novel_signatures:
        try:
            novel_signatures = pd.read_csv(args.novel_signatures_file, sep=',', header=0, index_col=0)
            _logger.info(f"\nNovel signature/s {[f'Manual{i}' for i in novel_signatures.index]} loaded: {novel_signatures.shape}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to load novel signatures file: {e}")

    if novel_signatures is not None:
        novel_signatures = novel_signatures.reindex(columns=COSMIC.columns, fill_value=0)
        novel_signatures.index = ["Manual" + str(i) for i in novel_signatures.index]
        COSMIC = pd.concat([COSMIC, novel_signatures], axis=0)

    if args.numpri == -1:
        num_runs = args.numruns
        current_primary = args.primary.copy()

        _logger.info(f'Recognition Phase is running..')
        recognition_results = run_stage(recognition_one_run, num_runs, args.numworkers, args, COSMIC, current_primary)
        list_primary = []
        for run_primary in recognition_results:
            list_primary.extend(run_primary)
            
        element_counts = Counter(list_primary)
        elements_with_weights = [(element, count / num_runs) for element, count in element_counts.items()]
        args.primary = [element for element, weight in elements_with_weights if weight >= args.thr3]
        _logger.debug(f'-> Weights: {elements_with_weights}')    
        _logger.info(f'-> Current Signature Set: {args.primary}')    
        
        _logger.info(f'Denovo Discovery Phase is running..')
        current_primary = args.primary.copy()
        while not stoppingCriteria: # as we have more novel signatures, it will continue.
            discovery_results = run_stage(discovery_one_run, num_runs, args.numworkers, args, COSMIC, current_primary)
            all_novel_profiles_list = []
            list_primary = []
            for run_primary, novel_profiles_one in discovery_results:
                list_primary.extend(run_primary)
                if not novel_profiles_one.empty:
                    all_novel_profiles_list.append(novel_profiles_one)
            
            if all_novel_profiles_list:
                all_novel_profiles_one = pd.concat(all_novel_profiles_list, ignore_index=True)
            else:
                all_novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
                
            element_counts = Counter(list_primary)
            elements_with_weights = [(element, count / num_runs) for element, count in element_counts.items()]
            args.primary = [element for element, weight in elements_with_weights if weight >= args.thr3]
            all_novel_profiles_one.to_csv(os.path.join(args.out,f'Collected_signature_profiles.txt'))
            
            if set(args.primary) == set(current_primary) and all_novel_profiles_one.shape[0] >= 2:
                _logger.debug(f'-> Weights: {elements_with_weights}')    
                _logger.info(f'-> Current signature set: {args.primary}')                               
                _logger.debug(f'{all_novel_profiles_one.shape[0]} novel signature profiles inferred.')
                X = all_novel_profiles_one[COSMIC.columns].values
                
                sim = cosine_similarity(X)  # pairwise cosine similarity matrix
                
                # adjacency: edge where similarity >= threshold (no self-edge)
                adj = (sim >= args.thr4).astype(int)
                np.fill_diagonal(adj, 0)
                
                # connected components (single-link at threshold)
                _, labels_cc = connected_components(csgraph=adj, directed=False, connection='weak')
                
                # gather components
                components = {}
                for idx, comp in enumerate(labels_cc):
                    components.setdefault(comp, []).append(idx)
                
                pruned_clusters = []
                for inds in components.values():
                    if len(inds) < 2:
                        continue
                    pruned = prune_component(inds, sim, args.thr4)
                    if len(pruned) >= 2:
                        pruned_clusters.append(pruned)
                
                # choose cluster
                chosen_indices = []
                if pruned_clusters:
                    chosen_indices = max(pruned_clusters, key=len)
                else:
                    # fallback to largest connected component
                    if components:
                        largest = max(components.values(), key=len)
                        if len(largest) >= 2:
                            chosen_indices = largest
                
                if chosen_indices:
                    chosen_labels = [all_novel_profiles_one.index[i] for i in chosen_indices]
                    _logger.debug(f"Cluster found with {len(chosen_indices)} profiles (threshold={args.thr4}): {chosen_labels}")
                    #cluster_medians = all_novel_profiles_one.loc[chosen_labels].median().to_frame().T
                    cluster_medians = (all_novel_profiles_one.iloc[chosen_indices][COSMIC.columns].median().to_frame().T)
                    # Normalize rows to sum to 1 (avoid division by zero)
                    row_sums = cluster_medians.sum(axis=1).replace(0, 1)
                    cluster_medians = cluster_medians.div(row_sums, axis=0)

                    cos = cosine_similarity(cluster_medians, COSMIC)[0]
                    best_idx = cos.argmax()
                    best_score = cos[best_idx]
                    best_name = COSMIC.index[best_idx]

                    if round(float(best_score), 2) > (args.thr7):
                        _logger.debug(f"Not a good novel candidate since it has cosine similarity {best_score:.2f} (> {args.thr7}) with {best_name}!")
                        chosen_indices = None
                                                    
                if chosen_indices and (len(chosen_indices) >= 2):
                    novel_count += 1
                    _logger.info(f'Checking confidence of Novel{novel_count}..')
                    args.out = current_args_out + f'/_withInferredNovel{novel_count}'
                    os.makedirs(args.out, exist_ok=True)
                    cluster_medians_path = os.path.join(args.out, f'Signature_profile_inferred.txt')
                    cluster_medians.index = ["Novel" + str(novel_count)]
                    trial_COSMIC = pd.concat([COSMIC, cluster_medians], axis=0)
                    trial_COSMIC.to_csv(os.path.join(args.out, f'COSMIC_now.txt'), index=True)
                    trial_primary = current_primary + ["Novel" + str(novel_count)]
                    
                    validation_results = run_stage(validation_one_run, num_runs, args.numworkers, args, trial_COSMIC, trial_primary)
                    list_primary = []
                    for run_primary in validation_results:
                        list_primary.extend(run_primary)
                    
                    element_counts = Counter(list_primary)
                    elements_with_weights = [(element, count / num_runs) for element, count in element_counts.items()]
                    validated_primary = [element for element, weight in elements_with_weights if weight >= args.thr3]
                    
                    if set(trial_primary) == set(validated_primary):
                        _logger.info(f"Novel{novel_count} accepted.")
                        cluster_medians.to_csv(cluster_medians_path)            
                        _logger.info(f"Novel{novel_count}'s signature profile saved at: {cluster_medians_path}")
                        stoppingCriteria = False
                        args.primary = trial_primary.copy()
                        current_primary = args.primary.copy()
                        combined_cluster_path = os.path.join(args.out, f'Signature_profiles_all.txt')
                       
                        if novel_signatures is not None:
                            combined_txt = pd.concat([novel_signatures, cluster_medians])
                        else:
                            combined_txt = cluster_medians
                        
                        combined_txt.to_csv(combined_cluster_path)
                        novel_signatures = combined_txt    
                        
                        args.novel_signatures_file = combined_cluster_path
                        COSMIC = trial_COSMIC
                    
                    else:
                        _logger.info(f"Novel candidate Novel{novel_count} rejected.")
                        stoppingCriteria = True
                        args.primary = current_primary.copy()
                        novel_count -= 1
            
                else:
                    _logger.info(f"Novel {novel_count + 1} is NOT accepted.")
                    stoppingCriteria = True
                    args.primary = current_primary
                    args.numpri = current_numpri
                
            else:
                _logger.info(f'No novel signature accepted!')
                args.primary = current_primary
                args.numpri = current_numpri
                stoppingCriteria = True  
            
    else:
        _logger.info(f"Inference with {args.numpri} novel signatures (no Redenovo recognition)..")
        
    if current_numpri == -1:
        args.numpri = 0
        
    args.out = current_args_out
    
    #initialize tensors
    data = resources.Resources(args, COSMIC)
    
    #initialize optimizer
    opt = model.Optimizer(data)
    
    #perform optimization
    best_loss_for_run = opt.optimize()
    
    _logger.debug(f'Final run with {args.primary}: (best loss = {best_loss_for_run:.2f})')
    _logger.info(f'Final run with {args.primary}')
     
    #store in resources
    opt.store()        
            
    #export to file
    export = exporter.Exporter(data)
    export.write_tables()
    export.write_tables_exposure()
    export.write_logs()
    _logger.info(f'DONE! Total time taken: {round(time.time()-start_time)} seconds')


def GetKnownSBSs(cosine_sim_df, d, threshold):
    known_SBSs = []
    weights = []
    
    for i in range(d.shape[0]):
        row_values = cosine_sim_df.iloc[i]
        top_3_columns = row_values.nlargest(3)
        
        for j in top_3_columns.index:
            val = cosine_sim_df.loc[cosine_sim_df.index[i], j]

            if hasattr(val, "__len__") and not isinstance(val, (float, int)):
                val = max(val)

            val = round(float(val), 2)

            if val >= threshold:
                known_SBSs.append(j)
                weights.append(val)

    new_df = pd.DataFrame({
        "SBS": known_SBSs,
        "Weight": weights
    })

    if not new_df.empty:
        new_df["SBS"] = new_df["SBS"].astype(str)
        new_df = new_df[~new_df["SBS"].str.startswith("ReDeNovo")]

    return new_df["SBS"].tolist(), new_df["Weight"].tolist()


# prune a component so each member has similarity >= threshold to ALL other members
def prune_component(indices, sim_matrix, thr):
    inds = list(indices)
    changed = True
    while changed and len(inds) > 1:
        changed = False
        to_remove = []
        for i in inds:
            others = [j for j in inds if j != i]
            if not others:
                continue
            min_sim = sim_matrix[i, others].min()
            if min_sim < thr:
                to_remove.append(i)
        if to_remove:
            for r in to_remove:
                inds.remove(r)
            changed = True
    return inds


def get_COSMIC(genome, WGS_or_WES, exclude, version):
    if WGS_or_WES == 'WES':
        file_path = os.path.join(tool_dir, 'COSMIC_catalogue', f'COSMIC_v{version}_SBS_GRCh{genome}_exome.txt')
    else:
        file_path = os.path.join(tool_dir, 'COSMIC_catalogue', f'COSMIC_v{version}_SBS_GRCh{genome}.txt')    
        
    sigDB = pd.read_table(file_path, header=0)
    COSMIC = sigDB.iloc[:, 1:].T.apply(pd.to_numeric)
    COSMIC.index = sigDB.columns[1:]
    COSMIC.columns = sigDB['Type']
    
    if len(exclude) > 0:
        missing_exclude = [sig for sig in exclude if sig not in COSMIC.index]
        if missing_exclude:
            raise ValueError(f"Excluded signature(s) not found in COSMIC: {missing_exclude}")
        COSMIC = COSMIC.drop(exclude)
    
    return COSMIC

            
def recognition_one_run(run_id):
    args = copy.copy(_WORKER_STATE["args"])
    cosmic_df = _WORKER_STATE["cosmic_df"]
    current_primary = list(_WORKER_STATE["current_primary"])
    args.primary = current_primary.copy()
    banHistory = Counter()
    banned_fixed = []
    
    def ordered_unique(items):
        return list(dict.fromkeys(items))
    
    def persistent_bans():
        return [s for s, c in banHistory.items() if c >= args.thr6]
    
    # Tracks the primary-set size at the last ban event.
    ban_primary_length = len(args.primary)
    initial_pass_completed = False
    i = 1
    SBS_df = pd.DataFrame(columns=['iter1'])

    while i <= args.numiters:
        args.numpri = i
        data = resources.Resources(args, cosmic_df)
        opt = model.Optimizer(data)
        opt.optimize()
        opt.store()

        normalized_data = data.A.div(data.A.sum(axis=1), axis=0)
        condition = normalized_data >= args.exposure_thr1
        percent = condition.mean()
        condition2 = data.A >= args.exposure_thr2
        percent2 = condition2.mean()

        if initial_pass_completed and i == 1:
            percentages = percent[~percent.index.str.startswith('ReDeNovo')]
            percentages2 = percent2[~percent2.index.str.startswith('ReDeNovo')]
            subdata_p = data.P['fixed'].copy()
            mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)

            if mask.sum() != subdata_p.shape[0]:
                failed_sigs = subdata_p.index[~mask].tolist()
                args.primary = subdata_p.index[mask].tolist()
        
                banHistory.update(failed_sigs)
        
                banned_fixed = ordered_unique(
                    banned_fixed + failed_sigs + persistent_bans()
                )
        
                ban_primary_length = len(args.primary)
        
                _logger.debug(
                    f"Run {run_id}: Fixed SBS set updated after exposure filtering: "
                    f"{args.primary}"
                )
        
                i = 1
                continue
        
            else:
                # clear temporary bans, but keep signatures that failed at least x times historically.
                if len(args.primary) > ban_primary_length:
                    banned_fixed = persistent_bans()
                    ban_primary_length = len(args.primary)

        initial_pass_completed = True
        percentages = percent[percent.index.str.startswith('ReDeNovo')]
        percentages2 = percent2[percent2.index.str.startswith('ReDeNovo')]
        subdata_p = data.P['inferred'].copy()
        mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
        subdata_p = subdata_p[mask]

        if subdata_p.shape[0] > 0:
            novels_only = subdata_p[subdata_p.index.str.startswith('ReDeNovo')]
            cosine_sim = cosine_similarity(novels_only, cosmic_df)
            cosine_sim_df = pd.DataFrame(cosine_sim, index=novels_only.index, columns=cosmic_df.index)
            new_SBSs, weights = GetKnownSBSs(cosine_sim_df, novels_only, args.thr2)
            
            candidate_pairs = [
                (sbs, weight)
                for sbs, weight in zip(new_SBSs, weights)
                if sbs not in args.primary and sbs not in banned_fixed
            ]
            
            _logger.debug(f"Run {run_id}: Current set: {args.primary}")
            _logger.debug(f"Run {run_id}: Suggested Reference Signature/s: {set(new_SBSs)}")
            _logger.debug(f"Run {run_id}: Suggested Reference Signature/s (excluding primary and banned {banned_fixed}): {set(sbs for sbs, _ in candidate_pairs)}")
                        
            if candidate_pairs:
                candidate_df = pd.DataFrame(candidate_pairs, columns=["SBS", "Weight"])
                candidate_df["Weight"] = pd.to_numeric(candidate_df["Weight"], errors="coerce")
                candidate_df = candidate_df.dropna()
                per_iter_weights = candidate_df.groupby("SBS")["Weight"].max()
            else:
                per_iter_weights = pd.Series(dtype=float)
            
            SBS_df = SBS_df.reindex(
                SBS_df.index.union(per_iter_weights.index),
                fill_value=0.0
            )
            
            column_name = f"iter{i}"
            SBS_df[column_name] = 0.0
            
            if not per_iter_weights.empty:
                SBS_df.loc[per_iter_weights.index, column_name] = per_iter_weights.astype(float)
        
            SBS_df_dedup = SBS_df[~SBS_df.index.duplicated(keep="first")]

            hit_counts = (SBS_df_dedup >= args.thr2).sum(axis=1)
            valid = SBS_df_dedup[hit_counts >= args.consno]
            
            if not valid.empty:
                avg_scores = valid.where(valid >= args.thr2).mean(axis=1)
                best_row = avg_scores.idxmax()
                args.primary.append(best_row)
                _logger.debug(f"Run {run_id} Added {best_row} to current signature set: {args.primary}")
            
                SBS_df = pd.DataFrame(columns=["iter1"])
                i = 1
                continue

        i += 1
    _logger.info(f"Run {run_id} inferred: {args.primary}")

    return args.primary

def discovery_one_run(run_id):
    args = copy.copy(_WORKER_STATE["args"])
    cosmic_df = _WORKER_STATE["cosmic_df"]
    current_primary = list(_WORKER_STATE["current_primary"])
    
    args.primary = current_primary.copy()
    args.numpri = 1

    data = resources.Resources(args, cosmic_df)
    opt = model.Optimizer(data)
    opt.optimize()
    opt.store()

    normalized_data = data.A.div(data.A.sum(axis=1), axis=0)
    condition = normalized_data >= args.exposure_thr1
    percent = condition.mean()
    condition2 = data.A >= args.exposure_thr2
    percent2 = condition2.mean()
    
    novel_profiles_one = pd.DataFrame(columns=cosmic_df.columns)
    percentages = percent[~percent.index.str.startswith('ReDeNovo')]
    percentages2 = percent2[~percent2.index.str.startswith('ReDeNovo')]
    subdata_p = data.P['fixed'].copy()
    mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
    
    if mask.sum() != subdata_p.shape[0]:
        args.primary = subdata_p[mask].index.tolist()

    percentages = percent[percent.index.str.startswith('ReDeNovo')]
    percentages2 = percent2[percent2.index.str.startswith('ReDeNovo')]
    subdata_p = data.P['inferred'].copy()
    mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
    subdata_p = subdata_p[mask]
    
    if subdata_p.shape[0] == 0:
        return args.primary, novel_profiles_one

    novels_only = subdata_p[subdata_p.index.str.startswith('ReDeNovo')]
        
    if novels_only.shape[0] == 0:
        return args.primary, novel_profiles_one

    cosine_sim = cosine_similarity(novels_only, cosmic_df)
    cosine_sim_df = pd.DataFrame(cosine_sim, index=novels_only.index, columns=cosmic_df.index)

    filtered_sim = cosine_sim_df.loc[:, ~cosine_sim_df.columns.str.startswith('ReDeNovo')]
    
    if filtered_sim.shape[1] == 0:
        novel_profiles_one = novels_only.reset_index(drop=True)

        return args.primary, novel_profiles_one


    top_labels = filtered_sim.idxmax(axis=1)
    vals = np.array([
        max(np.atleast_1d(cosine_sim_df.loc[row, label]))
        for row, label in zip(cosine_sim_df.index, top_labels)
    ])

    rounded_vals = np.round(vals.astype(float), 2)

    novel_mask = rounded_vals <= args.thr7

    novel_profiles_one = (
        novels_only.iloc[np.where(novel_mask)[0]]
        .reset_index(drop=True)
    )

    return args.primary, novel_profiles_one

def validation_one_run(run_id):
    args = copy.copy(_WORKER_STATE["args"])
    cosmic_df = _WORKER_STATE["cosmic_df"]
    trial_primary = list(_WORKER_STATE["current_primary"])
    
    args.primary = trial_primary.copy()
    args.numpri = 0
    
    data = resources.Resources(args, cosmic_df)
    opt = model.Optimizer(data)
    opt.optimize()
    opt.store()

    normalized_data = data.A.div(data.A.sum(axis=1), axis=0)
    condition = normalized_data >= args.exposure_thr1
    percent = condition.mean()
    condition2 = data.A >= args.exposure_thr2
    percent2 = condition2.mean()

    percentages = percent[~percent.index.str.startswith('ReDeNovo')]
    percentages2 = percent2[~percent2.index.str.startswith('ReDeNovo')]
    subdata_p = data.P['fixed'].copy()
    mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
    
    if mask.sum() != subdata_p.shape[0]:
        args.primary = subdata_p[mask].index.tolist()

    return args.primary


def get_available_cpus():
    """
    Return the number of CPUs actually available to this job.

    Priority:
    1. SLURM_CPUS_PER_TASK
    2. CPU affinity
    3. os.cpu_count()
    """
    slurm_cpus = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm_cpus is not None:
        try:
            return max(1, int(slurm_cpus))
        except ValueError:
            pass

    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:
        return max(1, os.cpu_count() or 1)


def get_max_workers(num_runs, requested_workers=None, reserve_cpus=0):
    """
    Determine how many worker processes to use.

    Priority:
    1. User-requested workers
    2. Available CPUs

    The final number of workers is capped by num_runs,
    since using more worker processes than runs provides no benefit.

    This function manages only process-level parallelism.
    It does not control native threads used internally by NumPy,
    OpenBLAS, TensorFlow, or other libraries.
    """
    available = max(1, get_available_cpus() - max(0, int(reserve_cpus)))

    if requested_workers is None or requested_workers < 1:
        workers = min(num_runs, available)
    else:
        workers = min(requested_workers, available)

    return max(1, min(workers, num_runs))


def _init_worker(args, cosmic_df, current_primary):
    """
    Runs once in each subprocess.
    Stores shared stage data in process-local globals.
    """
    _WORKER_STATE["args"] = copy.copy(args)
    _WORKER_STATE["cosmic_df"] = cosmic_df
    _WORKER_STATE["current_primary"] = list(current_primary)


def _run_in_worker(worker_fn, run_id):
    try:
        return worker_fn(run_id)
    except Exception:
        _logger.exception(f"Worker failed in run {run_id}")
        raise

def _get_mp_context():
    """
    Use spawn for normal script/module execution.
    Use fork for notebook/interactive testing on Linux.
    """
    main_mod = sys.modules.get("__main__")
    is_interactive = (
        "ipykernel" in sys.modules
        or main_mod is None
        or not hasattr(main_mod, "__file__")
    )

    if is_interactive and "fork" in mp.get_all_start_methods():
        _logger.warning(
            "Interactive session detected; using fork multiprocessing context. "
            "For production runs, use the installed command or python -m from a real .py module."
        )
        return mp.get_context("fork")

    return mp.get_context("spawn")

def run_stage(worker_fn, num_runs, requested_workers, args, cosmic_df, current_primary, reserve_cpus=0):
    max_workers = get_max_workers(num_runs, requested_workers, reserve_cpus)

    _logger.debug(
        f"Launching {max_workers} worker processes, "
        f"available_cpus={get_available_cpus()}"
    )

    # If only one worker process is requested, do NOT use ProcessPoolExecutor.
    if max_workers == 1:
        _init_worker(args, cosmic_df, current_primary)
        return [_run_in_worker(worker_fn, run) for run in range(num_runs)]

    results = [None] * num_runs
    ctx = _get_mp_context()

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=(args, cosmic_df, current_primary),
    ) as ex:
        futures = {
            ex.submit(_run_in_worker, worker_fn, run): run
            for run in range(num_runs)
        }

        for fut in as_completed(futures):
            run = futures[fut]
            results[run] = fut.result()

    return results


def print_config(args):
    print("\n========== ReDeNovo Configuration ==========")

    sections = {
        "Input": {
            "matrix": args.matrix,
            "genome": args.genome,
            "whole": args.whole,
            "cosmic_version": args.cosmic_version,
            "manual_cosmic": args.manual_cosmic,
            "manual_cosmic_file": args.manual_cosmic_file,
            "numruns": args.numruns,
            "numiters": args.numiters,
            "numworkers": args.numworkers,
            "delimiter": repr(args.delimiter),
            "verbosity": args.verbosity,
        },

        "Output": {
            "out": args.out,
        },

        "Primary Signatures": {
            "primary": args.primary,
            "exclude": args.exclude,
            "numpri": args.numpri,
        },
        
        "Novel Signatures": {
            "add_novel_signatures": args.add_novel_signatures,
            "novel_signatures_file": args.novel_signatures_file,
        },

        "Exposure Thresholds": {
            "exposure_thr1": args.exposure_thr1,
            "exposure_thr2": args.exposure_thr2,
            "exposure_thr3": args.exposure_thr3,
        },

        "Similarity / Presence Thresholds": {
            "thr1": args.thr1,
            "thr2": args.thr2,
            "thr3": args.thr3,
            "thr4": args.thr4,
            "thr5": args.thr5,
            "thr6": args.thr6,
            "thr7": args.thr7,
            "consno": args.consno,
        },

        "Optimizer": {
            "epochs": args.epochs,
            "stepsizes": args.stepsizes,
            "optimizers": args.optimizers,
        },

    }

    for section_name, values in sections.items():
        print(f"\n[{section_name}]")

        for key, value in values.items():
            print(f"{key:<30}: {value}")

    print("\n============================================\n")

def run():
    main(sys.argv[1:])
    
if __name__ == "__main__":
    run()
