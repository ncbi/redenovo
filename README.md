# ReDeNovo
Combined inference of known and novel mutational signatures with ReDeNovo

## Installation

Clone the GitHub repository and move into the project directory:

```bash
git clone https://github.com/ncbi/ReDeNovo.git
cd ReDeNovo

conda env create -f environment.yml
conda activate redenovo

pip install -e ./
redenovo -h
```

## Examples
Sample run:
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output
```
- specify genome, WGS/WES, and/or cosmic version (if COSMIC is used as catalogue)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -g 38 -w WGS --cosmic-version 3.4
```
- specify manual catalogue signature database (for modified COSMIC, subset of COSMIC, or any other catalogue)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output --manual-cosmic --manual-cosmic-file src/redenovo/data/COSMIC.txt
```
- parallelize (run with 10 cpus - see description of --numworkers to see how the number of cpus are determined)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output --numworkers 10 
```
- run with a given set of signatures (no denovo discovery phase, only refitting)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -N 0 -P SBS1 SBS2 SBS3 SBS5 SBS8 
```
- run with a given set of signatures (with given number of novel signatures, here 2)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -N 2 -P SBS1 SBS2 SBS3 SBS5 SBS8
```
- include a manual signature as catalogue and rerun the tool
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output --novel-signatures-file src/redenovo/manual_signature.txt --add-novel-signatures
```
- include a manual signature and rerun the tool with given set of signatures included (does not guarantee that given set will be kept, just force it to begin)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -P SBS1 SBS2 SBS3 SBS5 SBS8 --novel-signatures-file src/redenovo/manual_signature.txt --add-novel-signatures
```
- exclude signature/s from the given reference database
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -E SBS2 SBS40a  
```
- if input mutation count matrix file has row and column names
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --has-header-and-index True
```
- update number of run and iteration (for example: 20 internal runs and 10 iterations)
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output -n 20 -i 10
```
- be more selective to add a signature waiting for it to be selected repeatedly at least two times
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --consno 2
```

## Hyperparameters and Definitions

### User-defined Parameters

| Parameter | Description | Default |
|---------|-------------|---------|
| `--consno` | Minimum number of times a signature must be selected to be included in the inferred signature set | `1` |
| `-i`, `--numiters` | Maximum number of iterations allowed while attempting to add a new fixed signature (patience for novel signature detection) | `10` |
| `-n`, `--numruns` | Number of runs to repeat the analysis | `10` |
| `--thr1` | Minimum fraction of patients with exposure ≥ `thr1` required for a signature to be considered present | `0.1` |
| `--thr2` | Minimum cosine similarity to match a signature with a known COSMIC signature and include in the inferred set | `0.75` |
| `--thr3` | Minimum exposure weight for a signature to contribute to the final exposure profile (For example, a signature is accepted as inferred if it is present in ≥7 runs out of 10 when `--thr3` = 0.70 and `--numruns` = 10) | `0.70` |
| `--thr4` | Minimum cosine similarity required to connect two novel signature profiles during single-linkage clustering in the denovo discovery phase | `0.75` |
| `--thr5` | Minimum fraction of the cohort with nonzero exposure required for a signature to be considered present | `0.1` |
| `--thr6` | Number of times a signature can fail the support criteria within a run before being permanently banned for that run | `3` |
| `--thr7` | Minimum cosine similarity to consider a signature as known and exclude it from novel candidate detection | `0.75` |
| `--exposure-thr1` | Minimum patient-wise normalized exposure required for a signature to be considered present | `0.05` |
| `--exposure-thr2` | Minimum raw exposure required for a signature to be considered present | `1` |
| `-E`, `--exclude` | List of SBS signatures to exclude from COSMIC (e.g., `["SBS1", "SBS5"]`) | `[]` |
| `--numworkers` | Number of worker processes to use. If not provided or less than 1, uses min(numruns, available CPUs). Values greater than available CPUs are capped at the available CPU count. The final number of workers is always capped by numruns. | `1` |

### Database of Catalogue Signatures

| Parameter | Description | Default |
|---------|-------------|---------|
| `-g`, `--genome` | Genome version to use (`37` or `38`) | `38` |
| `-w`, `--whole` | Sequencing platform for the data and COSMIC catalogue (`WGS` or `WES`) | `WGS` |
| `--cosmic-version` | COSMIC version to use (`3.4`, `3.3`, `3.2`, `3.1`, `3`, `2`, or `1`) | `3.4` |
| `--manual-cosmic` | Whether to use a user-provided COSMIC reference file (`True` or `False`) | `False` |
| `--manual-cosmic-file` | Path to the file containing reference signatures (used only if `--manual-cosmic` is enabled) | `None` |

### Input File–Related Parameters

| Parameter | Description | Default |
|---------|-------------|---------|
| `-d`, `--delimiter` | Delimiter used to separate columns in input matrices; also used for output files | `\t` |
| `--has-header-and-index` | Whether the input file contains row names and column headers (`True` or `False`) | `False` |
| `--add-novel-signatures` | Whether to evaluate using novel signatures provided in an external file (`True` or `False`) | `False` |
| `--novel-signatures-file` | Path to the file containing novel signatures (used only if `--add-novel-signatures` is enabled) | `None` |

### Output Options

| Parameter | Description | Default |
|---------|-------------|---------|
| `-O`, `--out` | Path to the output directory. The folder will be created if it does not exist. Existing files will be overwritten. | Current directory |
