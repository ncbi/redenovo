# redenovo
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
