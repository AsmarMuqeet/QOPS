OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[5];
creg cregless[4];
h qregless[0];
h qregless[1];
h qregless[2];
h qregless[3];
x qregless[4];
cp(pi/16) qregless[4],qregless[0];
cp(pi/8) qregless[4],qregless[1];
id qregless[4];
cp(pi/4) qregless[4],qregless[2];
cp(pi/2) qregless[4],qregless[3];
swap qregless[0],qregless[3];
h qregless[0];
swap qregless[1],qregless[2];
cp(-pi/2) qregless[1],qregless[0];
h qregless[1];
cp(-pi/4) qregless[2],qregless[0];
cp(-pi/2) qregless[2],qregless[1];
h qregless[2];
cp(-pi/8) qregless[3],qregless[0];
cp(-pi/4) qregless[3],qregless[1];
cp(-pi/2) qregless[3],qregless[2];
h qregless[3];
barrier qregless[0],qregless[1],qregless[2],qregless[3],qregless[4];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
