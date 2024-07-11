OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[2];
creg cregless[2];
h qregless[1];
cp(pi) qregless[1],qregless[0];
barrier qregless[0],qregless[1];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
