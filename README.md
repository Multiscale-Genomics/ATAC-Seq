# ATAC-Seq

[![Documentation Status](https://readthedocs.org/projects/ATAC-Seq/badge/?version=latest)](http://ATAC-Seq.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Multiscale-Genomics/ATAC-Seq.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/ATAC-Seq) [![Code Health](https://landscape.io/github/Multiscale-Genomics/ATAC-Seq/master/landscape.svg?style=flat)](https://landscape.io/github/Multiscale-Genomics/ATAC-Seq/master)


This repository contains pipeline wrapper for the analysis of ATAC Seq data. ATAC Seq protocols are used for transposase accessible chromatin regions. 

# Requirements
- pyenv and pyenv-virtualenv
- Python 2.7.12
- Python Modules:
  - pylint
  - pytest
  - mg-tool-api
- Cutadapt
- Bowtie2
- Trim Galore
- R


Installation
------------

Directly from GitHub:

```
cd ${HOME}/code

git clone https://github.com/Multiscale-Genomics/ATAC-Seq.git

cd ATAC-Seq
```

Create the Python environment

```
pyenv-virtualenv 2.7.12 ATAC-Seq
pyenv activate ATAC-Seq
pip install -e .
pip install -r requirements.txt
```
