OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[5];
creg cregless[5];
ry(-5.48681214738515) qregless[0];
ry(6.17430358738484) qregless[1];
cz qregless[0],qregless[1];
ry(-2.48211471276918) qregless[2];
cz qregless[0],qregless[2];
cz qregless[1],qregless[2];
ry(-4.87073612147293) qregless[3];
cz qregless[0],qregless[3];
cz qregless[1],qregless[3];
cz qregless[2],qregless[3];
ry(-1.13225298444705) qregless[4];
cz qregless[0],qregless[4];
ry(1.48976006426089) qregless[0];
cz qregless[1],qregless[4];
ry(-5.70472180821892) qregless[1];
rzz(7*pi/4) qregless[0],qregless[1];
cz qregless[2],qregless[4];
ry(4.73160065476197) qregless[2];
cz qregless[0],qregless[2];
cz qregless[1],qregless[2];
cz qregless[3],qregless[4];
ry(-0.976675728274714) qregless[3];
cz qregless[0],qregless[3];
cz qregless[1],qregless[3];
cz qregless[2],qregless[3];
ry(3.94275633297612) qregless[4];
cz qregless[0],qregless[4];
ry(-0.930600976197059) qregless[0];
cz qregless[1],qregless[4];
ry(-4.09301167598876) qregless[1];
cz qregless[0],qregless[1];
cz qregless[2],qregless[4];
ry(-4.44929727040748) qregless[2];
cz qregless[0],qregless[2];
cz qregless[1],qregless[2];
cz qregless[3],qregless[4];
ry(-5.14922111308235) qregless[3];
cz qregless[0],qregless[3];
cz qregless[1],qregless[3];
cz qregless[2],qregless[3];
ry(5.36795345777456) qregless[4];
cz qregless[0],qregless[4];
ry(-4.94728127471714) qregless[0];
cz qregless[1],qregless[4];
ry(-0.166201936680721) qregless[1];
cz qregless[2],qregless[4];
ry(-4.11536885095106) qregless[2];
cz qregless[3],qregless[4];
ry(-0.893201142031406) qregless[3];
ry(-4.14563312095378) qregless[4];
barrier qregless[0],qregless[1],qregless[2],qregless[3],qregless[4];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
measure qregless[4] -> cregless[4];
