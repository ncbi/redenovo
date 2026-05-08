import logging
import sys, os
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
import time
from scipy.sparse.csgraph import connected_components
from redenovo import __version__
from redenovo import arguments
from redenovo import resources
from redenovo import optimizer as model
from redenovo import exporter

__author__ = "ReDeNovo"
__copyright__ = "ReDeNovo"
__license__ = "GPL-3.0-only"
_logger = logging.getLogger(__name__)

tool_dir = os.path.dirname(os.path.abspath(__file__))

def setup_logging(verbosity):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s: %(message)s"
    logging.basicConfig(
        level=logging.getLevelName(verbosity), stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S"
    )

# ---- CLI ----
# The functions defined in this section are wrappers around the main Python
# API allowing them to be called directly from the terminal as a CLI
# executable/script.

def main(args):
    print(f"ReDeNovo version: {__version__}")
    args = arguments.parse_args(args)
    os.makedirs(args.out, exist_ok=True)
    stoppingCriteria = False
    
    # setup the root logger
    setup_logging(args.verbosity)
    
    start_time = time.time()
    current_numpri = args.numpri
    current_args_out = args.out
    novel_count = 0
    min_link_threshold = 0.80
                    
    print(f'{args.matrix} {args.genome} {args.whole} init:{args.primary} exc:{args.exclude} consno:{args.consno} thr1:{args.thr1} thr2:{args.thr2} thr3:{args.thr3} thr4:{args.thr4} thr5:{args.thr5}\n')
    novel_signatures = None
    if args.manual_cosmic:
        if args.manual_cosmic_file is None:
            raise ValueError("You must provide --manual-cosmic-file when using --manual-cosmic.")
        
        try:
            COSMIC = pd.read_table(args.manual_cosmic_file, header=0)    
            COSMIC = COSMIC.iloc[:, 1:].T.apply(pd.to_numeric)
            print("Manual COSMIC loaded:", COSMIC.shape)

        except Exception as e:
            raise RuntimeError(f"Failed to load manual COSMIC file: {e}")

    else:
        COSMIC = get_COSMIC(args.genome, args.whole, args.exclude, args.cosmic_version)

    if args.add_novel_signatures:
        if args.novel_signatures_file is None:
            raise ValueError("You must provide --novel-signatures-file when using --add-novel-signatures.")
        
        try:
            novel_signatures = pd.read_csv(args.novel_signatures_file, sep=',', header=0, index_col=0)
            print(f"\nNovel signature/s {[f'Manual{i}' for i in novel_signatures.index]} loaded: {novel_signatures.shape}")
            
        except Exception as e:
            raise RuntimeError(f"Failed to load novel signatures file: {e}")

    if novel_signatures is not None:
        novel_signatures = novel_signatures.reindex(columns=COSMIC.columns, fill_value=0)
        novel_signatures.index = ["Manual" + str(i) for i in novel_signatures.index]
        COSMIC = pd.concat([COSMIC, novel_signatures], axis=0)
        
    # Recognition Phase
    if args.numpri == -1:
        num_runs = args.numruns
        num_iters = args.numiters
        cons_no = args.consno 
        list_primary = []
        all_novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
        current_primary = args.primary.copy()
        
        for run in range(num_runs):
            novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
            args.primary = current_primary.copy()
            _logger.info(f'Run {run} - Beginning signature set: {args.primary}')
            SBS_df = pd.DataFrame(columns=['iter1'])
            i = 1
            banned_fixed = []
            first_iter_flag = 1
            ban_length = 0
            
            while i <= num_iters:
                args.numpri = i
                column_name = 'iter' + str(args.numpri)
                
                data = resources.Resources(args, COSMIC)
                
                opt = model.Optimizer(data)
            
                best_loss_for_run = opt.optimize()
    
                opt.store()
                
                _logger.info(f'+{i}. iter: (best loss = {best_loss_for_run:.2f})')
                
                normalized_data = data.A.div(data.A.sum(axis=1), axis=0)
                condition = normalized_data >= args.exposure_thr1
                percent = condition.mean()
                condition2 = data.A >= args.exposure_thr2
                percent2 = condition2.mean()
                
                if first_iter_flag == 0 and i == 1:
                    percentages = percent[~percent.index.str.startswith('ReDeNovo')]
                    percentages2 = percent2[~percent2.index.str.startswith('ReDeNovo')]
                    subdata_p_inferred = data.P['fixed'].copy()
                    mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
                    
                    if sum(mask == True) != subdata_p_inferred.shape[0]:
                        args.primary = subdata_p_inferred[mask].index.tolist()
                        _logger.info(f'Updating set: {args.primary}')
                        i = 1
                        banned_fixed.extend(subdata_p_inferred[~mask].index.tolist())
                        ban_length = len(args.primary)
                        continue
                    
                    # if any signature is added
                    else:
                        if len(args.primary) > ban_length:
                            banned_fixed = []
                
                first_iter_flag = 0
                percentages = percent[percent.index.str.startswith('ReDeNovo')]
                percentages2 = percent2[percent2.index.str.startswith('ReDeNovo')]
                subdata_p_inferred = data.P['inferred'].copy()
                mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
                subdata_p_inferred = subdata_p_inferred[mask]
    
                if subdata_p_inferred.shape[0] > 0:
                    novels_only = subdata_p_inferred[subdata_p_inferred.index.str.startswith('ReDeNovo')]
                    cosine_sim = cosine_similarity(novels_only, COSMIC)
                    cosine_sim_df = pd.DataFrame(cosine_sim, index=novels_only.index, columns=COSMIC.index)
                    new_SBSs, weights = GetKnownSBSs(cosine_sim_df, novels_only, args.thr2)
                   
                    # Get novel signature profiles for further analysis
                    if i == 1: 
                        novel_rows = []

                        for ii in range(subdata_p_inferred.shape[0]):
                            row_values = cosine_sim_df.iloc[ii].drop(cosine_sim_df.columns[ii]) # Exclude self-comparison
                            row_values = row_values[~row_values.index.str.startswith('ReDeNovo')]
                            top_3_columns = row_values.nlargest(3)
                            j = top_3_columns.index[0]
                            val = cosine_sim_df.loc[cosine_sim_df.index[ii], j] 
                            if hasattr(val, "__len__") and not isinstance(val, (float, int)):
                                val = max(val)
                                
                            if round(float(val), 2) < args.thr4:
                                novel_rows.append(subdata_p_inferred.iloc[ii])
                                
                        novel_profiles_one = pd.DataFrame(novel_rows, columns=COSMIC.columns).reset_index(drop=True)
                    
                    new_SBSs = [item for item in new_SBSs if item not in args.primary]
                    new_SBSs = [item for item in new_SBSs if item not in banned_fixed]
                    _logger.info(f"Suggested catalogue signature/s: {set(new_SBSs)}")
                    
                    column_name = f'iter{i}'
                    weight_mapping = dict(zip(new_SBSs, weights))
                    SBS_df = SBS_df.reindex(SBS_df.index.union(new_SBSs), fill_value=0)
                    SBS_df[column_name] = SBS_df.index.map(weight_mapping).fillna(0).astype(float)

                    SBS_df_dedup = SBS_df[~SBS_df.index.duplicated(keep='first')]

                    counts = (SBS_df_dedup >= args.thr2).sum(axis=1)
                    
                    valid = SBS_df_dedup[counts >= cons_no]
                    
                    if not valid.empty:
                        avg_scores = valid.mean(axis=1)
                        best_row = avg_scores.idxmax()
                    else:
                        best_row = None
                        
                    # Add the best known SBS with the maximum average
                    if best_row is not None:
                        args.primary.append(best_row)                    
                        _logger.info(f"Adding {best_row} to signature set..")
                        _logger.info(f'Current signature set: {args.primary}')
                            
                        SBS_df = pd.DataFrame(columns=['iter1'])  # Start with one column for the first iter
                        i = 1
                    else:
                        i = i+1
        
                else:
                    i = i+1
            
            _logger.info(f'Run {run} is done. Inferred Signature Set: {args.primary}')

            list_primary.extend(args.primary)
            if novel_profiles_one.shape[0] > 0:
                all_novel_profiles_one = pd.concat([all_novel_profiles_one, novel_profiles_one], ignore_index=True)
            
            
        element_counts = Counter(list_primary)
        
        # List all elements with their weights
        elements_with_weights = [(element, count / num_runs) for element, count in element_counts.items()]
        args.primary = [element for element, weight in elements_with_weights if weight >= args.thr3]
       
        _logger.info(f'-> Weights: {elements_with_weights}')    
        _logger.info(f'-> Final Signature Set: {args.primary}')    
        _logger.info(f'Denovo Discovery Phase')
                           
        # Denovo Discovery Phase
        while not stoppingCriteria: # as we have more novel signatures, it will continue.
            stoppingCriteria = True        

            if all_novel_profiles_one.shape[0] > 0:
                
                if all_novel_profiles_one.shape[0] < 2:
                    print(f'{all_novel_profiles_one.shape[0]} novel signature profile (<2) inferred, no clustering!')
                    
                else: 
                    print(f'{all_novel_profiles_one.shape[0]} novel signature profiles inferred.')
                    X = all_novel_profiles_one.select_dtypes(include=['number']).values
                    
                    # --- single-link / min-link clustering with cosine similarity (threshold = 0.80) ---
                    sim = cosine_similarity(X)  # pairwise cosine similarity matrix
                    
                    # adjacency: edge where similarity >= threshold (no self-edge)
                    adj = (sim >= min_link_threshold).astype(int)
                    np.fill_diagonal(adj, 0)
                    
                    # connected components (single-link at threshold)
                    _, labels_cc = connected_components(csgraph=adj, directed=False, connection='weak')
                    
                    # gather components
                    components = {}
                    for idx, comp in enumerate(labels_cc):
                        components.setdefault(comp, []).append(idx)
                    
                    pruned_clusters = []
                    for comp_idx, inds in components.items():
                        if len(inds) < 2:
                            continue
                        pruned = prune_component(inds, sim, min_link_threshold)
                        if len(pruned) >= 2:
                            pruned_clusters.append(pruned)
                    
                    # choose cluster
                    chosen_indices = []
                    if pruned_clusters:
                        chosen_indices = max(pruned_clusters, key=len)
                    else:
                        if components:
                            largest = max(components.values(), key=len)
                            if len(largest) >= 2:
                                chosen_indices = largest
                    
                    # generate cluster_medians (element-wise median of chosen rows) and normalize
                    if chosen_indices:
                        chosen_labels = [all_novel_profiles_one.index[i] for i in chosen_indices]
                        print(f"Cluster found with {len(chosen_indices)} profiles (threshold={min_link_threshold}): {chosen_labels}")
                        cluster_medians = all_novel_profiles_one.loc[chosen_labels].median().to_frame().T
                        row_sums = cluster_medians.sum(axis=1).replace(0, 1)
                        cluster_medians = cluster_medians.div(row_sums, axis=0)
    
                        cos = cosine_similarity(cluster_medians, COSMIC)[0]
                        best_idx = cos.argmax()
                        best_score = cos[best_idx]
                        best_name = COSMIC.index[best_idx]
    
                        if round(float(best_score), 2) >= (args.thr4):
                            print(f"Not a good novel candidate since it has cosine similarity {best_score:.2f} (>= {args.thr4}) with {best_name}!")
                            chosen_indices = None
                        else:
                            # Assign simple Cluster labels for plotting: 0 for chosen cluster, -1 for others
                            all_novel_profiles_one['Cluster'] = -1
                            for idx in chosen_indices:
                                all_novel_profiles_one.iloc[idx, all_novel_profiles_one.columns.get_loc('Cluster')] = 0
                                                        
                    if chosen_indices and (len(chosen_indices) >= 2):
                        novel_count += 1
                        print(f'Checking confidence of Novel {novel_count}..')
                        args.out = current_args_out + f'/_withInferredNovel{novel_count}'
                        os.makedirs(args.out, exist_ok=True)
                        cluster_medians_path = os.path.join(args.out, f'Signature_profile_inferred.txt')
                        cluster_medians.index = ["Novel" + str(novel_count)]
                        COSMIC = pd.concat([COSMIC, cluster_medians], axis=0)
                        COSMIC.to_csv(os.path.join(args.out, f'COSMIC_now.txt'), index=True)
                        args.primary = args.primary + ["Novel" + str(novel_count)]
                        list_primary = []
                        all_novel_profiles_list = []
                        current_primary = args.primary.copy()
                        
                        for run in range(num_runs):
                            novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
                            args.primary = current_primary.copy()
                            
                            SBS_df = pd.DataFrame(columns=['iter1'])
                            i = 1
                            banned_fixed = []
                            first_iter_flag = 1
                            ban_length = 0
                            
                            args.numpri = i
                            column_name = 'iter' + str(args.numpri)
                            
                            data = resources.Resources(args, COSMIC)
                            
                            opt = model.Optimizer(data)
                        
                            best_loss_for_run = opt.optimize()
                
                            opt.store()
                            
                            normalized_data = data.A.div(data.A.sum(axis=1), axis=0)
                            condition = normalized_data >= args.exposure_thr1
                            percent = condition.mean()
                            condition2 = data.A >= args.exposure_thr2
                            percent2 = condition2.mean()
                            
                            if first_iter_flag == 0 and i == 1:
                                percentages = percent[~percent.index.str.startswith('ReDeNovo')]
                                percentages2 = percent2[~percent2.index.str.startswith('ReDeNovo')]
                                subdata_p_inferred = data.P['fixed'].copy()
                                mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
                                
                                if mask.sum() != subdata_p_inferred.shape[0]:
                                    args.primary = subdata_p_inferred[mask].index.tolist()
                                    i = 1
                                    banned_fixed.extend(subdata_p_inferred[~mask].index.tolist())
                                    ban_length = len(args.primary)
            
                                    continue
                                
                                # if any signature is added
                                else:
                                    if len(args.primary) > ban_length:
                                        banned_fixed = []
                                        
                            first_iter_flag = 0
                            percentages = percent[percent.index.str.startswith('ReDeNovo')]
                            percentages2 = percent2[percent2.index.str.startswith('ReDeNovo')]
                            subdata_p_inferred = data.P['inferred'].copy()
                            mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
                            subdata_p_inferred = subdata_p_inferred[mask]
                
                            if subdata_p_inferred.shape[0] > 0:
                                novels_only = subdata_p_inferred[subdata_p_inferred.index.str.startswith('ReDeNovo')]
                                cosine_sim = cosine_similarity(novels_only, COSMIC)
                                cosine_sim_df = pd.DataFrame(cosine_sim, index=novels_only.index, columns=COSMIC.index)
                                new_SBSs, weights = GetKnownSBSs(cosine_sim_df, novels_only, args.thr2)
                               
                                # Get novel signature profiles for further analysis
                                if i == 1: 
                                    novel_rows = []
            
                                    for ii in range(subdata_p_inferred.shape[0]):
                                        row_values = cosine_sim_df.iloc[ii].drop(cosine_sim_df.columns[ii]) # Exclude self-comparison
                                        row_values = row_values[~row_values.index.str.startswith('ReDeNovo')]
                                        top_3_columns = row_values.nlargest(3)
                                        j = top_3_columns.index[0]
                                        val = cosine_sim_df.loc[cosine_sim_df.index[ii], j] 
                                        if hasattr(val, "__len__") and not isinstance(val, (float, int)):
                                            val = max(val)
                                            
                                        if round(float(val), 2) < args.thr4:
                                            novel_rows.append(subdata_p_inferred.iloc[ii])
                                            
                                    novel_profiles_one = pd.DataFrame(novel_rows, columns=COSMIC.columns).reset_index(drop=True)
                        
                                new_SBSs = [item for item in new_SBSs if item not in args.primary]
                                new_SBSs = [item for item in new_SBSs if item not in banned_fixed]
                                
                                column_name = f'iter{i}'
                                weight_mapping = dict(zip(new_SBSs, weights))
                                SBS_df = SBS_df.reindex(SBS_df.index.union(new_SBSs), fill_value=0)
                                SBS_df[column_name] = SBS_df.index.map(weight_mapping).fillna(0).astype(float)
                                
                                SBS_df_dedup = SBS_df[~SBS_df.index.duplicated(keep='first')]
            
                                counts = (SBS_df_dedup >= args.thr2).sum(axis=1)
                                
                                valid = SBS_df_dedup[counts >= cons_no]
                                
                                if not valid.empty:
                                    avg_scores = valid.mean(axis=1)
                                    best_row = avg_scores.idxmax()
                                else:
                                    best_row = None

                                # Add the best known SBS with the maximum average
                                if best_row is not None:
                                    args.primary.append(best_row)                    
                                    
                            _logger.info(f'Run {run} is done. Inferred signature set: {args.primary}')
                
                            list_primary.extend(args.primary)
                            if not novel_profiles_one.empty:
                                all_novel_profiles_list.append(novel_profiles_one)

                        if all_novel_profiles_list:
                            all_novel_profiles_one = pd.concat(all_novel_profiles_list, ignore_index=True)
                        else:
                            all_novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
                            
                        # Count the frequency of each element in the list
                        element_counts = Counter(list_primary)
                        
                        # List all elements with their weights
                        elements_with_weights = [(element, count / num_runs) for element, count in element_counts.items()]                         
                        args.primary = [element for element, weight in elements_with_weights if weight >= args.thr3]
                        
                        all_novel_profiles_one.to_csv(os.path.join(args.out,f'Collected_signature_profiles.txt'))
                            
                        if set(args.primary) == set(current_primary + ["Novel" + str(novel_count)]):
                            stoppingCriteria = False
                            
                            _logger.info(f'-> Weights: {elements_with_weights}')    
                            _logger.info(f'-> Final signature set: {args.primary}')                               
                            cluster_medians.to_csv(cluster_medians_path)            
                            print(f"Novel {novel_count} is accepted and its signature profile saved at: {cluster_medians_path}")
                            combined_cluster_path = os.path.join(args.out, f'Signature_profiles_all.txt')
                           
                            if novel_signatures is not None:
                                combined_txt = pd.concat([novel_signatures, cluster_medians])#, ignore_index=True)
                                combined_txt.to_csv(combined_cluster_path)
                            else:
                                cluster_medians.to_csv(combined_cluster_path)
                                novel_signatures = cluster_medians 
                                
                            args.numpri = 1
                            args.novel_signatures_file = combined_cluster_path
                        else:
                            print(f"Novel {novel_count} is NOT accepted.")
                            stoppingCriteria = True
                            args.primary = current_primary
                            args.numpri = current_numpri
                    else:
                        print(f"No cluster with (>=2 profiles) satisfied the min-link constraint (>= {min_link_threshold}).")
     
            else:
                print(f'No novel signature inferred!')
                
    else:
        _logger.info(f"Run (wihout full approach) with {args.numpri} novel signatures to infer..")
        
    if current_numpri == -1:
        args.numpri = 0
        
    args.out = current_args_out
    
    data = resources.Resources(args, COSMIC)
    
    opt = model.Optimizer(data)
    
    best_loss_for_run = opt.optimize()
    
    _logger.info(f'Final run with {args.primary}: (best loss = {best_loss_for_run:.2f})')
     
    opt.store()        
            
    export = exporter.Exporter(data)
    export.write_tables()
    export.write_tables_exposure()
    export.write_logs()
    
    
    _logger.info(f'Total time: {time.time()-start_time} sec')

def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])

if __name__ == "__main__":
    run()


def GetKnownSBSs(cosine_sim_df, d, threshold):
    known_SBSs = []
    weights = []
    
    for i in range(d.shape[0]):
        row_values = cosine_sim_df.iloc[i]
        top_3_columns = row_values.nlargest(3)
        
        for j in top_3_columns.index[:3]:
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
    
    if len(exclude)>0:
        COSMIC = COSMIC.drop(exclude)
    
    return COSMIC
