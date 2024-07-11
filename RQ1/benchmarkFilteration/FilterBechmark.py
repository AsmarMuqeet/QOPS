#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from qiskit import QuantumCircuit, transpile, Aer
from gensim.matutils import hellinger
from forest.benchmarking.distance_measures import total_variation_distance
from scipy.spatial.distance import jensenshannon
from qiskit_ibm_runtime import Sampler, Session, QiskitRuntimeService, Options
import json
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import shutil
from tqdm import tqdm
service = QiskitRuntimeService(channel="ibm_quantum",token="69f16f44de2970998b5b3a8dc5dc8622760f117523728ab8201ad4f07013b2ae8829944e812b7f453bc0ff943ae345b52ed356255a30c686ab3f52a8de2ca894")
backend = service.get_backend('ibm_brisbane')
ideal_simulator = Aer.get_backend('qasm_simulator')


# In[ ]:


import multiprocessing
import subprocess
import os
import signal
import time

def task_function(program_file,backend):
    NOTDONE = True
    while NOTDONE:
        try:
            parts = program_file.split(".")[0].split("_")
            program_name = parts[0]
            qubits = parts[-1]
            qc = QuantumCircuit.from_qasm_file("benchmark/{}".format(program_file))
            qc.remove_final_measurements()
            qc.measure_all()
            qc_trans = transpile(qc,backend,optimization_level=0,seed_transpiler=42)
            NOTDONE = False
            if qc_trans.depth()<500:
                #print(program_file)
                shutil.copy("benchmark/{}".format(program_file),"filtered_benchmark/{}".format(program_file))
        except:
            continue

# # MQT Total 504 (32q)

# In[ ]:


if __name__ == "__main__":
    progress = tqdm(os.listdir("benchmark/"))
    for program_file in sorted(os.listdir("benchmark/"))[250::]:
        task_process = multiprocessing.Process(target=task_function,args=(program_file,backend))
        task_process.start() 
        task_process.join(timeout=180)
        if task_process.is_alive():
            print("Task did not complete within the timeout. Killing the process.")
            os.kill(task_process.pid, signal.SIGTERM)
        progress.update(1)
