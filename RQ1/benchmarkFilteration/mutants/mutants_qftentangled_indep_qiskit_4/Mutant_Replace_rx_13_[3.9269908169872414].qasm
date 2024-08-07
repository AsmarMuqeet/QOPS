OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[4];
creg cregless[4];
h qregless[3];
cx qregless[3],qregless[2];
cx qregless[2],qregless[1];
cx qregless[1],qregless[0];
h qregless[3];
cp(pi/2) qregless[3],qregless[2];
h qregless[2];
cp(pi/4) qregless[3],qregless[1];
cp(pi/2) qregless[2],qregless[1];
h qregless[1];
cp(pi/8) qregless[3],qregless[0];
cp(pi/4) qregless[2],qregless[0];
cp(pi/2) qregless[1],qregless[0];
rx(5*pi/4) qregless[0];
swap qregless[0],qregless[3];
swap qregless[1],qregless[2];
barrier qregless[0],qregless[1],qregless[2],qregless[3];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
