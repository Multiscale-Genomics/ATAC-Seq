# ATAC-Seq

[![Documentation Status](https://readthedocs.org/projects/mg-process-test/badge/?version=latest)](http://mg-process-test.readthedocs.io/en/latest/?badge=latest) [![Build Status](https://travis-ci.org/Multiscale-Genomics/ATAC-Seq.svg?branch=master)](https://travis-ci.org/Multiscale-Genomics/ATAC-Seq) [![Code Health](https://landscape.io/github/Multiscale-Genomics/ATAC-Seq/master/landscape.svg?style=flat)](https://landscape.io/github/Multiscale-Genomics/ATAC-Seq/master)


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

git clone https://github.com/Multiscale-Genomics/mg-process-test.git

cd mg-process-test
```

Create the Python environment

```
pyenv-virtualenv 2.7.12 mg-process-test
pyenv activate mg-process-test
pip install -e .
pip install -r requirements.txt
```
