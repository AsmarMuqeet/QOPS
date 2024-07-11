#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

# Importing standard Qiskit libraries
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import Aer
from qiskit_aer.primitives import Estimator as ideal_estimator
from qiskit_ibm_provider import IBMProvider
from qiskit_aer.noise.noise_model import NoiseModel
from qiskit.quantum_info import SparsePauliOp
from qiskit_ibm_runtime import Estimator,Session, QiskitRuntimeService, Options
import json
import os
from itertools import product
from tqdm.notebook import tqdm
from mitiq.zne.inference import LinearFactory,RichardsonFactory,PolyFactory,ExpFactory,AdaExpFactory
from mitiq import zne,cdr
from collections import defaultdict


provider = IBMProvider(token="")


# In[2]:


def get_noise_model(device_backend):
    
    noise = NoiseModel.from_backend(device_backend)
    noise_dict = noise.to_dict()
    for i,error in enumerate(noise_dict["errors"]):
        if "ecr" in error["operations"]:
            error["operations"] = ["cx" if x=="ecr" else x for x in error["operations"]]
    return noise_dict


def load_program_mutant_mix(program_name,mutant_number):
    program = "benchmarkFilteration/benchmark/"+program_name
    mutant_number = mutant_number.split("_")[-1]
    folder = program_name.replace(".qasm","")
    mutant_name = sorted(os.listdir(f"benchmarkFilteration/mutants/mutants_{folder}/"))[int(mutant_number)-1]
    mutant = f"benchmarkFilteration/mutants/mutants_{folder}/"+mutant_name
    return program,mutant

def hamiltonain_to_observable(program_name,mutant_number):
    program_name = program_name.replace(".qasm",".json")
    result_file = open(f"hamiltonian_result/{program_name}","r")
    result_data = json.load(result_file)
    hamiltonian = result_data[mutant_number]["faultyHamiltonian"]
    result_file.close()
    ham_keys = list(hamiltonian.keys())
    observable = hamiltonian[ham_keys[0]]*SparsePauliOp(ham_keys[0])
    for k in ham_keys[1::]:
        observable+=hamiltonian[k]*SparsePauliOp(k)
    return observable



def load_program_data(n_qubit,device_backend,ideal_simulator):
    job_data = []
    with open(f"RQ2_programs/rq2programs_{n_qubit}.json","r") as file:
        programs_and_mutants = json.load(file)

    for program_name,mutant_number,oracle in zip(programs_and_mutants["programs"],
                                                 programs_and_mutants["mutants"],
                                                 programs_and_mutants["ground"]):

        program_file,mutant_file = load_program_mutant_mix(program_name,mutant_number)
        
        
        fp = open(program_file,"r")
        qs = "".join([x for x in fp.readlines() if "creg" not in x and "measure" not in x])
        original_circuit = QuantumCircuit.from_qasm_str(qs)
        fp.close()
        
        fp = open(mutant_file,"r")
        qs = "".join([x for x in fp.readlines() if "creg" not in x and "measure" not in x])
        mutant_circuit = QuantumCircuit.from_qasm_str(qs)
        fp.close()
        
        mutant_circuit.remove_final_measurements()
        
        observable = hamiltonain_to_observable(program_name,mutant_number)
        
        ideal_expectation = ideal_simulator.run(transpile(original_circuit,basis_gates=['cx', 'id', 'rz', 'sx', 'x'],optimization_level=0,seed_transpiler=42),observable).result().values[0]
        
        transpiled_circuit_noise = transpile(mutant_circuit.copy(),backend=device_backend,basis_gates=['cx', 'id', 'rz', 'sx', 'x'],optimization_level=0,seed_transpiler=42)
        transpiled_circuit_ideal = transpile(mutant_circuit.copy(),basis_gates=['cx', 'id', 'rz', 'sx', 'x'],optimization_level=0,seed_transpiler=42)
        transpiled_observable_noise = observable.copy().apply_layout(transpiled_circuit_noise.layout)
        transpiled_observable_ideal = observable.copy()
        
        
        
        job_data.append({"ideal_expectation":ideal_expectation,
                         "mutant_circuit_noisy":transpiled_circuit_noise,
                         "mutant_circuit_ideal":transpiled_circuit_ideal,
                         "mutant_observable_noisy":transpiled_observable_noise,
                         "mutant_observable_ideal":transpiled_observable_ideal,
                         "program_name":program_name,
                         "mutant_number":mutant_number,
                         "ground":oracle})
        
    return job_data


class Mitiq_execution:
    def __init__(self,mitigation_name,observable,noise_observable):
        self.observable = observable.copy()
        self.noise_observable = noise_observable.copy()
        self.name = mitigation_name
        self.zne_configurations = {"Linear":LinearFactory(scale_factors=[1,3,5]),
                              "Richard":RichardsonFactory(scale_factors=[1,3,5]),
                              "Poly3":PolyFactory(scale_factors=[1,2,3,4,5],order=3),
                              "Poly4":PolyFactory(scale_factors=[1,2,3,4,5],order=4)
                                  }

        self.std_obs = None
    
    def execute_no_noise(self,circuit):
        if circuit.num_qubits>self.observable.num_qubits:
            return aer_estimator.run(circuit,self.noise_observable).result().values[0]
        else:
            return aer_estimator.run(circuit,self.observable).result().values[0]

    def execute_noise(self,circuit):
        if circuit.num_qubits>self.observable.num_qubits:
            return noise_sim.run(circuit,self.noise_observable).result().values[0]
        else:
            return noise_sim.run(circuit,self.observable).result().values[0]
    
    def execute_std_no_noise(self,circuit):
        if circuit.num_qubits<self.std_obs.num_qubits:
            obs = SparsePauliOp("Z"*circuit.num_qubits)
            return aer_estimator.run(circuit,obs).result().values[0]
        else:
            return aer_estimator.run(circuit,self.std_obs).result().values[0]

    def execute_std_noise(self,circuit):
        if circuit.num_qubits<self.std_obs.num_qubits:
            obs = SparsePauliOp("Z"*circuit.num_qubits)
            return noise_sim.run(circuit,obs).result().values[0]
        else:
            return noise_sim.run(circuit,self.std_obs).result().values[0]
        

    def get_standard_error(self,noise_sim,circuit,num_qubits,configuration):
        
        if self.name == "zne":

            total_circuit = circuit.compose(circuit.inverse()).copy()
            observable = SparsePauliOp("Z"*num_qubits).apply_layout(circuit.layout)
            self.std_obs = observable
            exps = []
            for x in range(10):
                exp = zne.execute_with_zne(total_circuit,self.execute_std_noise,factory=self.zne_configurations[configuration])
                exps.append(exp)

            npexps = np.array(exps)
            error = np.abs(1-npexps)
            std_error = np.std(error)/len(exps)
            return std_error
        elif self.name=="cdr":
            
            total_circuit = circuit.compose(circuit.inverse()).copy()
            observable = SparsePauliOp("Z"*num_qubits).apply_layout(circuit.layout)
            self.std_obs = observable
            exps = []
            for x in range(10):
                exp = cdr.execute_with_cdr(total_circuit.copy(),executor=self.execute_std_noise,simulator=self.execute_std_no_noise)
                exps.append(exp)

            npexps = np.array(exps)
            error = np.abs(1-npexps)
            std_error = np.std(error)/len(exps)
            return std_error

            
    def get_mitigation(self,name,circuit,runs,configuration):
        if name=="zne":
            return [zne.execute_with_zne(circuit,self.execute_noise,factory=self.zne_configurations[configuration]) for _ in range(runs)]
        elif name=="cdr":
            return [cdr.execute_with_cdr(circuit,executor=self.execute_noise,simulator=self.execute_no_noise) for _ in range(runs)]


# In[3]:


backend_name = "ibm_osaka"
device_backend = provider.get_backend(backend_name)

noise = NoiseModel.from_dict(get_noise_model(device_backend))
noise_sim = ideal_estimator(skip_transpilation=True,
                            #run_options={"seed":42,"shots":4000},
                            run_options={"shots":4000},
                            backend_options={"noise_model":noise,"coupling_map":device_backend.configuration().coupling_map},
                            approximation=False)

aer_estimator = ideal_estimator(skip_transpilation=True,
                                run_options={"seed":42,"shots":4000},
                                approximation=False)


# In[ ]:


if __name__=="__main__":
    for qubits in range(2,11):
        data = load_program_data(qubits,device_backend,aer_estimator)
        for index,program_data in enumerate(data):
            #if index>1:
            #    break
            cdr_mitiq_executor = Mitiq_execution("cdr",program_data["mutant_observable_ideal"],
                                                 program_data["mutant_observable_noisy"])

            std_err = cdr_mitiq_executor.get_standard_error(noise_sim,program_data["mutant_circuit_noisy"],qubits,None)
            mitigated_exps = cdr_mitiq_executor.get_mitigation("cdr",program_data["mutant_circuit_noisy"],10,None)
            result = {}
            if os.path.exists(r"./results/{}&cdr.json".format(backend_name)):
                file = open(r"./results/{}&cdr.json".format(backend_name),"r")
                result = json.load(file)
                file.close()

            file = open(r"./results/{}&cdr.json".format(backend_name),"w")
            temp = result.get(program_data["program_name"],{})
            temp[program_data["mutant_number"]] = {
                "ground":program_data["ground"],
                "ideal_exp":program_data["ideal_expectation"],
                "std_err":std_err,
                "mitigated_exps":mitigated_exps}
            result[program_data["program_name"]] = temp          
            json.dump(result,file)
            file.close()
        
            
        print(qubits)
        #break


# In[ ]:




