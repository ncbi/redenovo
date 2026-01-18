"""
This is a skeleton file that can serve as a starting point for a Python
console script. To run this script uncomment the following lines in the
``[options.entry_points]`` section in ``setup.cfg``::

    console_scripts =
         fibonacci = redenovo.skeleton:run

Then run ``pip install .`` (or ``pip install -e .`` for editable mode)
which will install the command ``fibonacci`` inside your current environment.

Besides console scripts, the header (i.e. until ``_logger``...) of this file can
also be used as template for Python modules.

Note:
    This skeleton file can be safely removed if not needed!

References:
    - https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
    - https://pip.pypa.io/en/stable/reference/pip_install
"""

# import configargparse
import logging
import sys, os
# import pkg_resources
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
import time

from redenovo import __version__
from redenovo import arguments
from redenovo import resources
from redenovo import optimizer as model
from redenovo import exporter


__author__ = "ReDeNovo"
__copyright__ = "ReDeNovo"
__license__ = "GPL-3.0-only"

_logger = logging.getLogger(__name__)

tool_dir = os.path.dirname(os.path.abspath(__file__)) # '/gpfs/gsfs12/users/kesimogluz2/ReNovo/test_data/'

# ---- Python API ----
# The functions defined in this section can be imported by users in their
# Python scripts/interactive interpreter, e.g. via
# `from redenovo.skeleton import check_positive_int`,
# when using this Python module as a library.

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
    begin_time = time.time()
    args = arguments.parse_args(args)
    print(args)
    #output_dir = Path(args.out)
    #output_dir.mkdir(parents=True, exist_ok=True)
    
    final_text = f'{args.matrix} {args.genome} {args.whole} init:{args.primary} exc:{args.exclude} consno:{args.consno} thr1:{args.thr1} thr2:{args.thr2} thr3:{args.thr3} thr4:{args.thr4} thr5:{args.thr5}'
    if args.manual_cosmic:
        if args.manual_cosmic_file is None:
            raise ValueError("You must provide --manual-cosmic-file when using --manual-cosmic.")
        
        # Read the file
        try:
            COSMIC = pd.read_table(args.manual_cosmic_file, header=0)    
            COSMIC = COSMIC.iloc[:, 1:].T.apply(pd.to_numeric)
            COSMIC.index = sigDB.columns[1:]
            COSMIC.columns = sigDB['Type']
            if len(args.exclude)>0:
                COSMIC = COSMIC.drop(args.exclude)
        
            print("Manual COSMIC loaded:", COSMIC.shape)

        except Exception as e:
            raise RuntimeError(f"Failed to load manual COSMIC file: {e}")

    else:
        COSMIC = get_COSMIC(args.genome, args.whole, args.exclude, args.cosmic_version)
    
    # setup the root logger
    setup_logging(args.verbosity)
    start_time = time.time()

    novel_signatures = None
    if args.add_novel_signatures:
        if args.novel_signatures_file is None:
            raise ValueError("You must provide --novel-signatures-file when using --add-novel-signatures.")
        
        # Read the file
        try:
            novel_signatures = pd.read_csv(args.novel_signatures_file, sep=',', header=0, index_col=0)
            print("Novel signatures loaded:", novel_signatures.shape)
            final_text = final_text + '\n' + '"Novel signatures loaded:", {novel_signatures.shape}'
        
        except Exception as e:
            raise RuntimeError(f"Failed to load novel signatures file: {e}")

    if novel_signatures is not None:
        # Align columns
        novel_signatures = novel_signatures.reindex(columns=COSMIC.columns, fill_value=0)
        novel_signatures.index = ["Novel" + str(i) for i in novel_signatures.index]
        # Append the new signature/s
        print(novel_signatures.index)
        COSMIC = pd.concat([COSMIC, novel_signatures], axis=0)
        
    
    if args.numpri == -1:
        num_runs = args.numruns
        num_iters = args.numiters
        cons_no = args.consno 
        list_primary = []
        
        all_novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
        all_novel_profiles_two = pd.DataFrame(columns=COSMIC.columns)
        current_primary = args.primary.copy()
        
        for run in range(num_runs):
            novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
            novel_profiles_two = pd.DataFrame(columns=COSMIC.columns)
        
            args.primary = current_primary.copy()
            
            _logger.info(f'Results for {args.matrix} -> run {run}: Genome: {args.genome}  Seq: {args.whole}')
            _logger.info(f'Beginning Signature Set: {args.primary}')
    
            # dataframe to keep results for each iteration of the current run
            SBS_df = pd.DataFrame(columns=['iter1'])
            i = 1
            banned_fixed = []
            start_iter_time = time.time()
            first_iter_flag = 1
            ban_length = 0
            
            while i <= num_iters:
                args.numpri = i
                column_name = 'iter' + str(args.numpri)
                
                #initialize tensors
                data = resources.Resources(args, COSMIC)
                
                #initialize optimizer
                opt = model.Optimizer(data)
            
                #perform optimization
                best_loss_for_run = opt.optimize()
    
                #store in resources
                opt.store()
                
                _logger.info(f'+{i}. iter: (best loss = {best_loss_for_run:.2f}):')
                
                # Patient-wise normalized data (sum to 1)
                normalized_data = data.A.div(data.A.sum(axis=1), axis=0)
                #normalized_data[normalized_data.index.str.startswith('ReDeNovo')]
                # Check if each value in normalized DataFrame is >= 0.05
                condition = normalized_data >= args.exposure_thr1
                # Calculate the percentage of rows in each column that are >= 1 (nonzero)
                percent = condition.mean()
                condition2 = data.A >= args.exposure_thr2
                percent2 = condition2.mean()
                percent.name = 'perc_samples_exceeding_thr_raw'
                percent2.name = 'perc_samples_exceeding_thr_norm'
                print(pd.merge(percent, percent2, left_index=True, right_index=True))
   
                if first_iter_flag == 0 and i == 1:
                    percentages = percent[~percent.index.str.startswith('ReDeNovo')]
                    percentages2 = percent2[~percent2.index.str.startswith('ReDeNovo')]
                    subdata_p_inferred = data.P['fixed'].copy()
                    mask = (percentages >= args.thr1) & (percentages2 >= args.thr5)
                    
                    if sum(mask == True) != subdata_p_inferred.shape[0]:
                        args.primary = subdata_p_inferred[mask].index.tolist()
                        _logger.info(f'**Fixed SBS set is updated for {i}.iter: {args.primary}')
                        i = 1
                        banned_fixed.extend(subdata_p_inferred[-mask].index.tolist())
                        ban_length = len(args.primary)

                        print(f'***********************banned_fixed: {banned_fixed}, ban_length = {ban_length}')
                        continue
                    
                    # if any signature is added
                    else:
                        if len(args.primary) > ban_length:
                            banned_fixed = []
                            print(f'*********************** -> banned_fixed: {banned_fixed}, ban_length = {ban_length}')
                        
                
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
                    new_SBSs, weights, ids = GetKnownSBSs(cosine_sim_df, novels_only, args.thr2)
                   
                    # Get novel signature profiles for further analysis
                    if i == 1 or i == 2: 
                        # just keep novel signatures from last time
                        if i == 1:
                            novel_profiles_one = pd.DataFrame(columns=COSMIC.columns)
                        if i == 2:
                            novel_profiles_two = pd.DataFrame(columns=COSMIC.columns)
        
                        
                        for ii in range(subdata_p_inferred.shape[0]):
                            row_values = cosine_sim_df.iloc[ii].drop(cosine_sim_df.columns[ii]) # Exclude self-comparison
                            row_values = row_values[~row_values.index.str.startswith('ReDeNovo')]
                            top_3_columns = row_values.nlargest(3)
                            j = top_3_columns.index[0]
                            val = cosine_sim_df.loc[cosine_sim_df.index[ii], j] 
                            
                            if round(val, 2) < args.thr4:
                                row = subdata_p_inferred.iloc[ii]
                                if i == 1:
                                    novel_profiles_one = pd.concat([novel_profiles_one, row.to_frame().T], ignore_index=True)
                                elif i == 2:
                                    novel_profiles_two = pd.concat([novel_profiles_two, row.to_frame().T], ignore_index=True)
                                
                                _logger.info(f"Current novel option for i={i}: {subdata_p_inferred.index[ii]}")
                                
                                      
                        
                    
                     
                    new_SBSs = [item for item in new_SBSs if item not in args.primary]
                    _logger.info(f"Suggested Reference Signature/s: {set(new_SBSs)}")
                    new_SBSs = [item for item in new_SBSs if item not in banned_fixed]
                    _logger.info(f"Suggested Reference Signature/s (excluding banned: {banned_fixed}): {set(new_SBSs)}")
                    
                    column_name = f'iter{i}'
                    weight_mapping = dict(zip(new_SBSs, weights))
                    # Check if each new_SBS is already in the DataFrame
                    missing_rows = [sbs for sbs in new_SBSs if sbs not in SBS_df.index]
                    if missing_rows:
                        missing_df = pd.DataFrame(0, index=missing_rows, columns=[column_name])
                        SBS_df = pd.concat([SBS_df, missing_df])
                    
                    SBS_df[column_name] = SBS_df.index.to_series().apply(lambda x: weight_mapping.get(x, 0)).astype(float)
                    print(SBS_df)
            
                    # Iterate through the DataFrame to find rows with any {cons_no} values >= threshold
                    best_row = None
                    max_avg = -1  
                    
                    for row in SBS_df.index:
                        row_values = SBS_df.loc[~SBS_df.index.duplicated(keep='first')].loc[row, :]  # row_values = SBS_df.loc[row, :]
                    
                        # Filter values above the threshold
                        # already above threshold, so removed: row_values = row_values[row_values >= args.thr2]
                        
                        # If there are at least x values above the threshold, calculate the mean
                        if len(row_values) >= cons_no:
                            avg = row_values.mean()  # Calculate the mean of those values
        
                            # Update the best row if this average is higher than the previous max
                            if avg > max_avg:
                                max_avg = avg
                                best_row = row
                                    
                                
                    # Add the best known SBS with the maximum average
                    if best_row:
                        args.primary.append(best_row)                    
                        _logger.info(f"Adding {best_row} to P_fixed - which has the maximum average ({max_avg:.2f}) [for any {cons_no} values >= {args.thr2}]")
                        _logger.info(f'Current Signature Set: {args.primary}')
                            
                        SBS_df = pd.DataFrame(columns=['iter1'])  # Start with one column for the first iter
                        i = 1
                    else:
                        i = i+1
        
                        
                else:
                    i = i+1
    
                current_iter_time = time.time()
                elapsed_iter_time = current_iter_time - start_iter_time
                start_iter_time = current_iter_time
                _logger.info(f'It took {elapsed_iter_time:.2f} sec for {i-1}. iter of {run}. run.')
    
            
            _logger.info(f'Run {run} is done. Inferred Signature Set: {args.primary}')

            list_primary = list_primary + args.primary
            if novel_profiles_one.shape[0] > 0:
                all_novel_profiles_one = pd.concat([all_novel_profiles_one, novel_profiles_one], ignore_index=True)
            
            if novel_profiles_two.shape[0] > 0:
                all_novel_profiles_two = pd.concat([all_novel_profiles_two, novel_profiles_two], ignore_index=True)
    
            current_time = time.time()
            elapsed_time = current_time - start_time
            start_time = current_time
            _logger.info(f'It took {elapsed_time:.2f} sec for {run}. run.')
       
        # Count the frequency of each element in the list
        _logger.info(f'\n All final selected for {num_runs} runs: {list_primary}')
        
        element_counts = Counter(list_primary)
        
        # List all elements with their weights
        elements_with_weights = [(element, count / num_runs) for element, count in element_counts.items()]
         
        args.primary = [element for element, weight in elements_with_weights if weight >= args.thr3]
       
        _logger.info(f'-> Weights: {elements_with_weights}')    
        _logger.info(f'-> Final Signature Set: {args.primary}')    

        all_weights_are_one = all(weight == 1.0 for element, weight in elements_with_weights)
        
        if all_weights_are_one:
            final_text = final_text + '\n' + f'{all_weights_are_one}! sel: {args.primary} with {len(args.primary)} SBSs'
        else:
            final_text = final_text + '\n' + f'Weights: {elements_with_weights} sel: {args.primary} with {len(args.primary)} SBSs'
    
        if all_novel_profiles_one.shape[0] > 0:
            all_novel_profiles_one.to_csv(os.path.join(args.out,'Novel_signature_profiles_one.txt'))#,#_'+ str(self.args.numpri) + '.txt'),
        
        if all_novel_profiles_two.shape[0] > 0:
            all_novel_profiles_two.to_csv(os.path.join(args.out,'Novel_signature_profiles_two.txt'))#,#_'+ str(self.args.numpri) + '.txt'),
    
        args.numpri = 0


    else:
        _logger.info(f"Last phase wihout Redenovo's fixed set adjustment with {args.numpri} novel signatures to infer..")
        final_text = final_text + '\n' + f"Last phase wihout Redenovo's fixed set adjustment with {args.numpri} novel signatures to infer.."
    
    
    print(args.primary)

    #initialize tensors
    data = resources.Resources(args, COSMIC)
    
    #initialize optimizer
    opt = model.Optimizer(data)
    
    #perform optimization
    best_loss_for_run = opt.optimize()
    
    _logger.info(f'Final run: (best loss = {best_loss_for_run:.2f}):')
    final_text = final_text + '\n' + f'best loss = {best_loss_for_run:.2f}\n\n'
    
    #store in resources
    opt.store()        
            
    #export to file
    export = exporter.Exporter(data)
    export.write_tables()
    export.write_tables_exposure()
    export.write_logs()

    
    current_time = time.time()
    elapsed_time = current_time - start_time
    start_time = current_time
    _logger.info(f'{elapsed_time} for final run.')
    
    _logger.info(f'Total time: {time.time()-begin_time} sec')
    print(final_text)
    
    with open(os.path.join(os.getcwd(),f'_redenovo_results_summary.txt'), "a") as file:
        file.write(final_text)



def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])

if __name__ == "__main__":
    # ^  This is a guard statement that will prevent the following code from
    #    being executed in the case someone imports this file instead of
    #    executing it as a script.
    #    https://docs.python.org/3/library/__main__.html

    # After installing your project with pip, users can also run your Python
    # modules as scripts via the ``-m`` flag, as defined in PEP 338::
    #
    #     python -m redenovo.skeleton 42
    #
    run()



def GetKnownSBSs(cosine_sim_df, d, threshold):
    known_SBSs = []
    weights = []
    ids = []
    text = ''
    for i in range(d.shape[0]):
        # Exclude self-comparison
        row_values = cosine_sim_df.iloc[i]
        # OLD block
        #row_values = cosine_sim_df.iloc[i].drop(cosine_sim_df.columns[i])
        #row_values = row_values[~row_values.index.str.startswith('ReDeNovo')]
        
        top_3_columns = row_values.nlargest(3)
        j = top_3_columns.index[0]
        val = cosine_sim_df.loc[cosine_sim_df.index[i], j] 
        
        if round(val, 2) >= threshold:
            known_SBSs.append(j)  # Append the most similar row index/name
            weights.append(round(val, 2))
            ids.append(i)
        
        text = text + f"{d.index[i]} => most similar to {j} with cos. sim. {round(val, 2)}, "
        
        j = top_3_columns.index[1]
        val = cosine_sim_df.loc[cosine_sim_df.index[i], j] 
        text = text + f"second: {j} with {round(val, 2)}, "
        
        if round(val, 2) >= threshold:
            known_SBSs.append(j)  # Append the second most similar row index/name
            weights.append(round(val, 2))
            ids.append(i)
            
        j = top_3_columns.index[2]
        val = cosine_sim_df.loc[cosine_sim_df.index[i], j]
        text = text + f"third: {j} with {round(val, 2)} \n"
        
        if round(val, 2) >= threshold:
            known_SBSs.append(j)  # Append the third most similar row index/name
            weights.append(round(val, 2))
            ids.append(i)
        
    print(text)    
    new_df = pd.DataFrame({
        'SBS': known_SBSs,
        'Weight': weights,
        'ID': ids
    })

    new_df['SBS'] = new_df['SBS'].astype(str)
    
    # Filter out rows if out of COSMIC (that is, exclude redenovo's novel signatures)
    new_df = new_df[~new_df['SBS'].str.startswith('ReDeNovo')]
    
    return(new_df['SBS'].tolist(), new_df['Weight'].tolist(), new_df['ID'].tolist())


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