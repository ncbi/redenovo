# ReDeNovo
Hybrid detection of known and novel mutational signatures with ReDeNovo

## Installation

Clone the GitHub repository and move into the project directory:

```bash
git clone https://github.com/ncbi/ReDeNovo.git
cd redenovo

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
- specify cosmic version, genome, and WGS/WES
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -g 38 -w WGS --cosmic-version 3.4
```
- specify manual reference signature database
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output --manual-cosmic --manual-cosmic-file src/redenovo/data/COSMIC.txt
```
- run with a given set of signatures (no novel signature inference)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -N 0 -P SBS1 SBS2 SBS3 SBS5 SBS8 SBS13 SBS17b SBS18 
```
- run with a given set of signatures (with given number of novel signature inference, here 2)
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -N 2 -P SBS1 SBS2 SBS3 SBS5 SBS8 SBS13 SBS17b SBS18 
```
- include the inferred signature and rerun the tool
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output --novel-signatures-file src/redenovo/output/Cluster_signature_profiles_one_k1.txt --add-novel-signatures
```
- include the inferred signature and rerun the tool with given set of signatures included
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -P SBS1 SBS2 SBS3 SBS5 SBS8 SBS13 SBS17b SBS18  --novel-signatures-file src/redenovo/output/Cluster_signature_profiles_one_k1.txt --add-novel-signatures
```
- exclude signature/s from the given reference database
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output -E SBS2  
```
- if input mutation count matrix file has row and column names
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --has-header-and-index True
```
- update number of run and iteration (for example: 20 runs and 10 iterations)
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --has-header-and-index True -n 20 -i 10
```
- be more selective to add a signature waiting for it to be selected repeatedly at least two times
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --has-header-and-index True --consno 2
```

## Hyperparameters and Definitions

### User-defined Parameters

| Parameter | Description | Default |
|---------|-------------|---------|
| `--consno` | Minimum number of times a signature must be selected to be included in the inferred signature set | `1` |
| `-i`, `--numiters` | Maximum number of iterations allowed while attempting to add a new fixed signature (patience for novel signature detection) | `10` |
| `-n`, `--numruns` | Number of runs to repeat the analysis | `10` |
| `--thr1` | Minimum fraction of patients with exposure ≥ `thr1` required for a signature to be considered present | `0.1` |
| `--thr2` | Minimum cosine similarity required to match a signature with a known COSMIC signature and include it in the inferred set | `0.70` |
| `--thr3` | Minimum exposure weight required for a signature to contribute to the final exposure profile | `0.70` (present in ≥7 runs out of 10) |
| `--thr4` | Minimum cosine similarity to consider a signature as known and exclude it from novel candidate detection | `0.80` |
| `--thr5` | Minimum fraction of the cohort with nonzero exposure required for a signature to be considered present | `0.1` |
| `--exposure-thr1` | Minimum patient-wise normalized exposure required for a signature to be considered present | `0.05` |
| `--exposure-thr2` | Minimum raw exposure required for a signature to be considered present | `1` |
| `-E`, `--exclude` | List of SBS signatures to exclude from COSMIC (e.g., `["SBS1", "SBS5"]`) | `[]` |

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
| `--check-novel` | Whether to check novel signatures (`True` or `False`) | `False` |
| `--novel-signatures-file` | Path to the file containing novel signatures (used only if `--add-novel-signatures` is enabled) | `None` |

### Output Options

| Parameter | Description | Default |
|---------|-------------|---------|
| `-O`, `--out` | Path to the output directory. The folder will be created if it does not exist. Existing files will be overwritten. | Current directory |
