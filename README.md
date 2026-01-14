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
- including the inferred signature and rerun the tool
```bash
redenovo -M src/redenovo/data/M.txt -O src/redenovo/output --novel-signatures-file src/redenovo/output/Cluster_signature_profiles_one_k1.txt --add-novel-signatures
```
- including the inferred signature and rerun the tool with given set of signatures are included
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
- update # run and iteration to 20 runs and 10 iterations
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --has-header-and-index True -n 20 -i 10
```
- to be more selective to add a signature waiting for it to be selected repeatedly at least two times
```bash
redenovo -M src/redenovo/data/redenovo_M.txt -O src/redenovo/output --has-header-and-index True --consno 2
```
