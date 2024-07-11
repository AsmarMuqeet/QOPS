# QOPS: Quantum Program Testing Through Commuting Pauli Strings on IBM’s Quantum Computers
This repository contains the code necessary to reproduce the results of QOPS.

# Installation
RQ1 is only tested on Windows 11 and RQ2 only supports Linux

### Dependencies

Anaconda Python distribution is required [here](https://www.anaconda.com/products/distribution):

Steps:

    1. Clone repository
    2. cd QOPS
    3  conda env create -f RQ1.yml
    4  conda env create -f RQ2.yml
    4. conda activate (qoin_rq1 or qoin_rq2) 
	

# Results:
### Downalod Raw result data for each RQ [here](https://drive.google.com/drive/folders/1focAA1vQ-N31YaMgxG9nGzzgBgbN_QiA?usp=sharing) and extract each RQ result to respective folders. After extraction there should be "hamiltonian_result" folder in each RQ folder.
RQ1.ipynb contains processed results for RQ1 and RQ2.ipynb contains results for RQ2

# Evaluate Research Questions:
### To Re-Evaluate RQ1 on new data
    python EvaluateQOPS.py
    
> **_NOTE:_** Evaluating RQs on new data might take several days, Depending on the system specifications.

### To Re-Evaluate RQ2 on new data
#### Run python files named as ErrorMitigation_Computer.py. For example (cdr_brisbane.py)

> **_NOTE:_** All files for RQ2 require an IBM API token for running on real computers.
