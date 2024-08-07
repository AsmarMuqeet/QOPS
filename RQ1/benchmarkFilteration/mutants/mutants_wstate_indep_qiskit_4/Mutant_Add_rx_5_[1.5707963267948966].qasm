OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[4];
creg cregless[4];
ry(-pi/4) qregless[0];
ry(-0.9553166181245093) qregless[1];
ry(-pi/3) qregless[2];
x qregless[3];
cz qregless[3],qregless[2];
rx(pi/2) qregless[2];
ry(pi/3) qregless[2];
cz qregless[2],qregless[1];
ry(0.9553166181245093) qregless[1];
cz qregless[1],qregless[0];
ry(pi/4) qregless[0];
cx qregless[2],qregless[3];
cx qregless[1],qregless[2];
cx qregless[0],qregless[1];
barrier qregless[0],qregless[1],qregless[2],qregless[3];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
