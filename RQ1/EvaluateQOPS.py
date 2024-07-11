#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer
from qiskit.quantum_info import hellinger_distance, SparsePauliOp
import json
import numpy as np
import pandas as pd
import os
import shutil
import itertools
from tqdm.auto import tqdm
#from psfam_update.pauli_organizer import *
from psfam.pauli_organizer import *
from scipy.stats import chisquare

ideal_simulator = Aer.get_backend('qasm_simulator')


# In[ ]:


def load_program(name):
    
    qc = QuantumCircuit.from_qasm_file("benchmarkFilteration/benchmark/{}".format(name))
    qc.remove_final_measurements()
    if len(qc.clbits)>0:
        for i in range(len(qc.clbits)):
            qc.measure(i,i)
    else:
        qc.measure_all()
        
    return qc.copy()


def load_mutants(name):
    
    qc_list = []
    for mutant in sorted(os.listdir("benchmarkFilteration/mutants/mutants_{}/".format(name.replace(".qasm","")))):
        qc = QuantumCircuit.from_qasm_file("benchmarkFilteration/mutants/mutants_{}/{}".format(name.replace(".qasm",""),mutant))
        qc.remove_final_measurements()
        if len(qc.clbits)>0:
            for i in range(len(qc.clbits)):
                qc.measure(i,i)
        else:
            qc.measure_all()
        
        qc_list.append(qc.copy())
        
    return qc_list

def load_zero_input(name):
    file = open("benchmarkFilteration/zero_input/{}".format(name.replace(".qasm",".json")),"r")
    inp = json.load(file)
    file.close()
    return inp


def get_ps(qc,simulator,prob=True):
    
    counts = simulator.run(transpile(qc,basis_gates=simulator.operation_names,
                                     optimization_level=0,seed_transpiler=42),
                           shots=4000,seed_simulator=42).result().get_counts()
    if isinstance(counts,list):
        PS = []
        for count in counts:
            temp = dict(sorted(count.items(),key=lambda x: int(x[0],2)))
            if prob:
                for k in temp.keys():
                    temp[k] = temp[k]/4000
            PS.append(temp)
        return PS
    else:
        ps = dict(sorted(counts.items(),key=lambda x: int(x[0],2)))
        if prob:
            for k in ps.keys():
                ps[k] = ps[k]/4000
        return ps
    

def get_random_hamoltonian(fam):
    pauli_dict = {}
    commuting_pauli = fam.to_string()
    totalpauli = random.randint(2,len(commuting_pauli)-1)
    selected_pauli = random.sample(commuting_pauli,totalpauli)
    for s in selected_pauli:
        pauli_dict[s] = np.random.random()
    return pauli_dict

    
def get_n_hamoltonian_coeff(qubits,famindex):
    PO   = PauliOrganizer(qubits)
    fam = PO.get_families()[::-1][famindex]
    Hamiltonian_decomp = get_random_hamoltonian(fam)
    for s in Hamiltonian_decomp.keys():
        PO.input_pauli_decomp(s,Hamiltonian_decomp[s])

    f = PO.calc_coefficients()
    fam1 = f[::-1][famindex]
    return {"pauli":fam1.to_string(),"coeff":fam1.get_coefficients(),"decomp":Hamiltonian_decomp,"family":fam1}
        



def get_exp_from_family(counts,cfs):
    measurement = 0
    for k in counts.keys():
        measurement = measurement+ counts[k]*(cfs[int(k,2)]/4000)
        
    return measurement


def qucat_oracle(ps,observed):
    A = set(ps.keys())
    B = set(observed.keys())
    if len(A.difference(B))!=0 or len(B.difference(A))!=0:
        return 0
    else:
        _,pvalue = chisquare(list(observed.values()),list(ps.values()))
        if pvalue<0.01:
            return 0
        else:
            return 1
        

    

def evaluate(simulator):
    
    programs_name_list = sorted(os.listdir("benchmarkFilteration/benchmark/"))
    
    for program_name in tqdm(programs_name_list):
        
        
        if int(program_name.split("_")[-1].replace(".qasm",""))>11:
            continue
        
        
        mutants = load_mutants(program_name)
        result = {}
        
        mutant_number = 0
        for mutant in tqdm(mutants,leave=False):
            
            qubits = mutant.num_qubits
            
            PO  = PauliOrganizer(qubits)
            family  = PO.get_families()[::-1]
            
            for famindex,fam in enumerate(family):
                
                mutant_counts = get_ps(mutant,simulator,prob=False)
                err = 0
                faulty_hamiltonian = {}
                ham_result = 1
                for hcount in range(1,1001):
                
                    hamiltonian = get_n_hamoltonian_coeff(qubits,famindex)
                    
                    qc = load_program(program_name)
                    hamiltonian["family"].apply_to_circuit(qc)
                    ps = get_ps(qc,simulator,prob=False)
                    hamiltonian["family"].apply_to_circuit(mutant)
                    mutant_counts = get_ps(mutant,simulator,prob=False)
                    
                    true_exp = get_exp_from_family(ps,hamiltonian["coeff"])
                    observed_exp = get_exp_from_family(mutant_counts,hamiltonian["coeff"])
                    err = np.round(np.abs(true_exp-observed_exp),3)
                    qucat_result = qucat_oracle(ps,mutant_counts)
                    if err>0.01:
                        ham_result = 0
                        faulty_hamiltonian = hamiltonian["decomp"]
                        break
                    elif qucat_result==1:
                        break

                mutant_number+=1
            
                result["mutant_{}".format(mutant_number)] = {"err":err,"hcount":hcount,
                                                             "faultyHamiltonian":faulty_hamiltonian,
                                                             "qucat_oracle":qucat_result,"ham_oracle":ham_result}
        
        file = open("hamiltonian_result/{}".format(program_name.replace(".qasm",".json")),"w")
        json.dump(result,file)
        file.close()
        
        
if __name__ == '__main__':
    evaluate(ideal_simulator)