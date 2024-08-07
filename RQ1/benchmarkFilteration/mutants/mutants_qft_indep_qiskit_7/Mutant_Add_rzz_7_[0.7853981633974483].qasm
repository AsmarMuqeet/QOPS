OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[7];
creg cregless[14];
h qregless[6];
cp(pi/2) qregless[6],qregless[5];
h qregless[5];
cp(pi/4) qregless[6],qregless[4];
cp(pi/2) qregless[5],qregless[4];
h qregless[4];
cp(pi/8) qregless[6],qregless[3];
rzz(pi/4) qregless[5],qregless[3];
cp(pi/4) qregless[5],qregless[3];
cp(pi/2) qregless[4],qregless[3];
h qregless[3];
cp(pi/16) qregless[6],qregless[2];
cp(pi/8) qregless[5],qregless[2];
cp(pi/4) qregless[4],qregless[2];
cp(pi/2) qregless[3],qregless[2];
h qregless[2];
cp(pi/32) qregless[6],qregless[1];
cp(pi/16) qregless[5],qregless[1];
cp(pi/8) qregless[4],qregless[1];
cp(pi/4) qregless[3],qregless[1];
cp(pi/2) qregless[2],qregless[1];
h qregless[1];
cp(pi/64) qregless[6],qregless[0];
cp(pi/32) qregless[5],qregless[0];
cp(pi/16) qregless[4],qregless[0];
cp(pi/8) qregless[3],qregless[0];
cp(pi/4) qregless[2],qregless[0];
cp(pi/2) qregless[1],qregless[0];
h qregless[0];
swap qregless[0],qregless[6];
swap qregless[1],qregless[5];
swap qregless[2],qregless[4];
barrier qregless[0],qregless[1],qregless[2],qregless[3],qregless[4],qregless[5],qregless[6];
measure qregless[0] -> cregless[7];
measure qregless[1] -> cregless[8];
measure qregless[2] -> cregless[9];
measure qregless[3] -> cregless[10];
measure qregless[4] -> cregless[11];
measure qregless[5] -> cregless[12];
measure qregless[6] -> cregless[13];
