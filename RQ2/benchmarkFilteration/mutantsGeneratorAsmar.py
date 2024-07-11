import tkinter as tk
from tkinter import filedialog
import os
import random
import pandas as pd

from qiskit import QuantumCircuit


import warnings

from qiskit.circuit import CircuitInstruction, Instruction

warnings.simplefilter(action='ignore', category=FutureWarning)


AllGates = ("x", "h", "p", "t", "s", "z", "y", "id", "rx", "ry", "rz", "sx", "swap", "rzz", "rxx", "cx", "cz", "ccx", "cswap", "cp") #All gates that are implemented
OneQubit = ("x", "h", "p", "t", "s", "z", "y", "id", "rx", "ry", "rz", "sx") #Gates that only affect one qubit
TwoQubit = ("swap", "rzz", "rxx", "cx", "cz", "cp") #Gates that affect two qubit
MoreThanTwoQubit = ("ccx", "cswap") #Gates that affect more than two qubit
OneParam = ("p", "rx", "ry", "rz", "rzz", "rxx", "cp") #Gates that affect the phase and needs to specify a phase
TwoParam = () #Gates that affect the phase and needs to specify a phase
ThreeParam = ("u") #Gates that affect the phase and needs to specify a phase
phaseValues = [0.0, 0.7853981633974483, 1.5707963267948966, 2.356194490192345, 3.141592653589793, 3.9269908169872414, 4.71238898038469, 5.497787143782138]


def getPositionGate(qc):
    column_names = ['Position', 'Gate', 'Params', 'Qubits']
    df = pd.DataFrame(columns=column_names)
    x = 0
    for instruction in qc.data:
        if (instruction[0].name != 'measure') and (instruction[0].name != 'barrier'):
            new_line = {'Position': x, 'Gate': instruction[0].name, 'Params': instruction[0].params, 'Qubits': instruction.qubits}
            new_df = pd.DataFrame.from_dict(new_line, orient='index').T
            df = pd.concat([df, new_df], ignore_index=True)
            x = x + 1
    return df


def getNewQubits(n_qubits, origin_qubits):
    if n_qubits == len(origin_qubits):
        new_qubits = origin_qubits
    else:
        if n_qubits < len(origin_qubits):
            new_qubits = origin_qubits[(len(origin_qubits)-n_qubits-1):-1]

    return new_qubits


def getAddInstructions(position_gate, gates):
    column_names = ['Name', 'Operator', 'Position', 'Gate', 'New_gate', 'Params', 'New_params', 'Qubits', 'New_qubits']
    df = pd.DataFrame(columns=column_names)

    for index, row in position_gate.iterrows():
        for gate in gates:
            new_qubits = 0
            new_params = []
            if gate in OneParam:
                new_param = random.choice(phaseValues)
                new_params.append(new_param)

            if gate in TwoQubit:
                if len(row.Qubits) >= 2:
                    new_qubits = getNewQubits(2, row.Qubits)
            elif gate in MoreThanTwoQubit:
                if len(row.Qubits) >= 3:
                    new_qubits = getNewQubits(3, row.Qubits)
            else:
                if len(row.Qubits) >= 1:
                    new_qubits = getNewQubits(1, row.Qubits)

            if new_qubits != 0:
                name = 'Mutant_Add_' + str(gate) + '_' + str(row.Position) + '_' + str(new_params) + '.qasm'
                new_instruction = {'Name': name, 'Operator': 'Add', 'Position': row.Position, 'Gate': 'Gap', 'New_gate': gate,
                                   'Params': 'Gap', 'New_params': new_params, 'Qubits': 'Gap', 'New_qubits': new_qubits}
                new_df = pd.DataFrame.from_dict(new_instruction, orient='index').T
                df = pd.concat([df, new_df], ignore_index=True)

    return df


def getRemoveInstructions(position_gate):
    column_names = ['Name', 'Operator', 'Position', 'Gate', 'New_gate', 'Params', 'New_params', 'Qubits', 'New_qubits']
    df = pd.DataFrame(columns=column_names)
    for index, row in position_gate.iterrows():
        name = 'Mutant_' + 'Remove' + '_' + str(row.Gate) + '_' + str(row.Position) + '.qasm'
        new_instruction = {'Name': name, 'Operator': 'Remove', 'Position': row.Position, 'Gate': row.Gate, 'New_gate': 'Gap', 'Params': row.Params, 'New_params': row.Params, 'Qubits': row.Qubits, 'New_qubits': 'Gap'}
        new_df = pd.DataFrame.from_dict(new_instruction, orient='index').T
        df = pd.concat([df, new_df], ignore_index=True)

    return df


def getReplaceInstructions(position_gate, gates):
    column_names = ['Name', 'Operator', 'Position', 'Gate', 'New_gate', 'Params', 'New_params', 'Qubits', 'New_qubits']
    df = pd.DataFrame(columns=column_names)
    for index, row in position_gate.iterrows():
        for new_gate in gates:
            new_qubits = 0
            if new_gate != row.Gate:
                new_params = []
                if new_gate in OneParam:
                    new_param = random.choice(phaseValues)
                    new_params.append(new_param)

                if (new_gate in TwoQubit and len(row.Qubits) == 2) or (new_gate in MoreThanTwoQubit and len(row.Qubits) == 3) or (new_gate in OneQubit) and (len(row.Qubits) == 1):
                    new_qubits = row.Qubits

                if new_qubits != 0:
                    name = 'Mutant_' + 'Replace' + '_' + str(new_gate) + '_' + str(row.Position) + '_' + str(new_params) + '.qasm'
                    new_instruction = {'Name': name, 'Operator': 'Replace', 'Position': row.Position, 'Gate': row.Gate, 'New_gate': new_gate, 'Params': row.Params, 'New_params': new_params, 'Qubits': row.Qubits, 'New_qubits': new_qubits}
                    new_df = pd.DataFrame.from_dict(new_instruction, orient='index').T
                    df = pd.concat([df,new_df], ignore_index=True)

    return df

def getPosibleInstructions(position_gate_df, qubits):
    if len(qubits) > 2:
        gates = AllGates
    elif len(qubits) > 1:
        gates = OneQubit + TwoQubit
    else:
        gates = OneQubit

    column_names = ['Name', 'Operator', 'Position', 'Gate', 'New_gate', 'Params', 'New_params', 'Qubits', 'New_qubits']
    df = pd.DataFrame(columns=column_names)

    df_add = getAddInstructions(position_gate_df, gates)
    df = pd.concat([df, df_add], ignore_index=True)

    df_remove = getRemoveInstructions(position_gate_df)
    df = pd.concat([df, df_remove], ignore_index=True)

    df_replace = getReplaceInstructions(position_gate_df, gates)
    df = pd.concat([df, df_replace], ignore_index=True)

    return df


def createGate(instruction):

    new_name = instruction.New_gate
    new_params = instruction.New_params
    new_qubits = instruction.New_qubits
    new_num_qubits = len(instruction.New_qubits)

    new_gate = CircuitInstruction(operation=Instruction(name=new_name, num_qubits=new_num_qubits, num_clbits=0, params=new_params), qubits=new_qubits)

    return new_gate


def createMutant(instruction, qc):
    mutant = QuantumCircuit(qc.qubits, qc.clbits)
    mutated = False
    for x, gate in enumerate(qc.data):
        if x == instruction.Position:
            if instruction.Operator == 'Add':
                new_gate = createGate(instruction)
                mutant.append(new_gate)
                mutant.append(gate)
            elif instruction.Operator == 'Replace':
                new_gate = createGate(instruction)
                mutant.append(new_gate)
            mutated = True
        else:
            mutant.append(gate)
    if mutated == False: #NEED TO CHECK THE QUBIT NUMBER WHERE TO ADD AT THE END OF CIRCUIT
        new_gate = createGate(instruction)
        mutant.append(new_gate)  # NEED TO CHECK THE PARAMETERS, QUBITS AND VALUES
    return mutant

def createMutants(mutants_path, df_instructions, selection, qc):
    if selection['All']:
        for index, row in df_instructions.iterrows():
            mutant = createMutant(row, qc)
            qasm_code = mutant.qasm()

            # Save QASM code to a file
            file_path = mutants_path + "/" + str(row.Name)
            with open(file_path, "w") as qasm_file:
                qasm_file.write(qasm_code)


def start(origin_file, save_path, selection):

    origin_filename = origin_file.split('/')[-1]
    origin_filename = origin_filename.split('.')[0]

    qc = QuantumCircuit.from_qasm_file(origin_file)

    position_gates_df = getPositionGate(qc)
    df_instructions = getPosibleInstructions(position_gates_df, qc.qubits)
    #df_instructions.to_csv(save_path + '/instructions.csv', index=False)

    mutants_path = save_path + '/mutants_' + origin_filename
    os.mkdir(mutants_path)
    createMutants(mutants_path, df_instructions, selection, qc)

    pass


if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    # Prompt the user to select a .qasm file using a file dialog
    files = filedialog.askopenfilenames(title="Select the origin .qasm file", filetypes=[("QASM files", "*.qasm")])
    folder = filedialog.askdirectory(title="Select a dict to save the generated mutants")

    characteristics = {'All': True}

    for file in files:
        start(file, folder, characteristics)
