OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[5];
creg cregless[5];
h qregless[0];
h qregless[1];
cz qregless[0],qregless[1];
h qregless[2];
cz qregless[1],qregless[2];
h qregless[3];
cz qregless[0],qregless[3];
h qregless[4];
cz qregless[3],qregless[4];
barrier qregless[0],qregless[1],qregless[2],qregless[3],qregless[4];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
measure qregless[4] -> cregless[4];
