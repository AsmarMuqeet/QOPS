// Benchmark was created by MQT Bench on 2023-06-29
// For more information about MQT Bench, please visit https://www.cda.cit.tum.de/mqtbench/
// MQT Bench version: v1.0.0
// Qiskit version: {'qiskit-terra': '0.24.1', 'qiskit-aer': '0.12.0', 'qiskit-ignis': None, 'qiskit-ibmq-provider': '0.20.2', 'qiskit': '0.43.1', 'qiskit-nature': '0.6.2', 'qiskit-finance': '0.3.4', 'qiskit-optimization': '0.5.0', 'qiskit-machine-learning': '0.6.1'}

OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg meas[3];
h q[0];
h q[1];
rzz(6.84212273097045) q[0],q[1];
h q[2];
rzz(6.84212273097045) q[0],q[2];
rx(-0.42073950596148635) q[0];
rzz(6.84212273097045) q[1],q[2];
rx(-0.42073950596148635) q[1];
rzz(-6.153399127068727) q[0],q[1];
rx(-0.42073950596148635) q[2];
rzz(-6.153399127068727) q[0],q[2];
rx(-6.497223513925254) q[0];
rzz(-6.153399127068727) q[1],q[2];
rx(-6.497223513925254) q[1];
rx(-6.497223513925254) q[2];
barrier q[0],q[1],q[2];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
