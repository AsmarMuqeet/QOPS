OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[10];
creg cregless[10];
u3(1.4368347742816616,1.1646500873100205,-pi) qregless[0];
u3(0.13038834331032634,-0.2928382424047804,0) qregless[1];
cx qregless[0],qregless[1];
u3(2.3018560275705338,-3.1167849646094092,-pi) qregless[2];
cx qregless[0],qregless[2];
cx qregless[1],qregless[2];
u3(1.5783117544539513,0.07660625016677525,-pi) qregless[3];
cx qregless[0],qregless[3];
cx qregless[1],qregless[3];
cx qregless[2],qregless[3];
u3(3.1322119352256292,-1.1773372206208812,0) qregless[4];
cx qregless[0],qregless[4];
cx qregless[1],qregless[4];
cx qregless[2],qregless[4];
cx qregless[3],qregless[4];
u3(1.4124389803026798,-2.434570523812673,0) qregless[5];
cx qregless[0],qregless[5];
cx qregless[1],qregless[5];
cx qregless[2],qregless[5];
cx qregless[3],qregless[5];
cx qregless[4],qregless[5];
u3(1.2444656817555668,-1.7482629013133648,0) qregless[6];
cx qregless[0],qregless[6];
cx qregless[1],qregless[6];
cx qregless[2],qregless[6];
cx qregless[3],qregless[6];
cx qregless[4],qregless[6];
cx qregless[5],qregless[6];
u3(1.5046299106322627,-1.3076812305427237,-pi) qregless[7];
cx qregless[0],qregless[7];
cx qregless[1],qregless[7];
cx qregless[2],qregless[7];
cx qregless[3],qregless[7];
cx qregless[4],qregless[7];
cx qregless[5],qregless[7];
cx qregless[6],qregless[7];
u3(1.0625547235745711,-0.5166404252966217,0) qregless[8];
cx qregless[0],qregless[8];
cx qregless[1],qregless[8];
cx qregless[2],qregless[8];
cx qregless[3],qregless[8];
cx qregless[4],qregless[8];
cx qregless[5],qregless[8];
cx qregless[6],qregless[8];
cx qregless[7],qregless[8];
u3(0.5550554224571163,-1.7933732440688743,0) qregless[9];
cx qregless[0],qregless[9];
u3(2.8742785055981956,1.917773905749895,-pi) qregless[0];
cx qregless[1],qregless[9];
u3(0.8932807542109367,-3.005579583727833,0) qregless[1];
cx qregless[0],qregless[1];
cx qregless[2],qregless[9];
u3(2.345769178126651,-0.5739760098973861,0) qregless[2];
cx qregless[0],qregless[2];
cx qregless[1],qregless[2];
cx qregless[3],qregless[9];
u3(2.0474788819188667,-1.135773149735491,-pi) qregless[3];
cx qregless[0],qregless[3];
rx(3*pi/2) qregless[1];
cx qregless[1],qregless[3];
cx qregless[2],qregless[3];
cx qregless[4],qregless[9];
u3(2.7761197097590844,0.5683728542359914,0) qregless[4];
cx qregless[0],qregless[4];
cx qregless[1],qregless[4];
cx qregless[2],qregless[4];
cx qregless[3],qregless[4];
cx qregless[5],qregless[9];
u3(2.72699034602209,1.8893541777246625,0) qregless[5];
cx qregless[0],qregless[5];
cx qregless[1],qregless[5];
cx qregless[2],qregless[5];
cx qregless[3],qregless[5];
cx qregless[4],qregless[5];
cx qregless[6],qregless[9];
u3(2.401640904800445,-2.4254077858804965,-pi) qregless[6];
cx qregless[0],qregless[6];
cx qregless[1],qregless[6];
cx qregless[2],qregless[6];
cx qregless[3],qregless[6];
cx qregless[4],qregless[6];
cx qregless[5],qregless[6];
cx qregless[7],qregless[9];
u3(3.059042641009883,2.0651656802006935,-pi) qregless[7];
cx qregless[0],qregless[7];
cx qregless[1],qregless[7];
cx qregless[2],qregless[7];
cx qregless[3],qregless[7];
cx qregless[4],qregless[7];
cx qregless[5],qregless[7];
cx qregless[6],qregless[7];
cx qregless[8],qregless[9];
u3(2.1966192898367836,-2.846934388642458,-pi) qregless[8];
cx qregless[0],qregless[8];
cx qregless[1],qregless[8];
cx qregless[2],qregless[8];
cx qregless[3],qregless[8];
cx qregless[4],qregless[8];
cx qregless[5],qregless[8];
cx qregless[6],qregless[8];
cx qregless[7],qregless[8];
u3(2.5067461861055573,0.7934855547557511,-pi) qregless[9];
cx qregless[0],qregless[9];
u3(2.8426000178928454,-0.6752586753862846,-pi) qregless[0];
cx qregless[1],qregless[9];
u3(1.1354532936221056,-2.554363801359381,-pi) qregless[1];
cx qregless[0],qregless[1];
cx qregless[2],qregless[9];
u3(1.2500242582094414,-1.1240263022165689,0) qregless[2];
cx qregless[0],qregless[2];
cx qregless[1],qregless[2];
cx qregless[3],qregless[9];
u3(0.8994360763247731,-2.1918765046211157,-pi) qregless[3];
cx qregless[0],qregless[3];
cx qregless[1],qregless[3];
cx qregless[2],qregless[3];
cx qregless[4],qregless[9];
u3(2.2094986973106154,2.413462260298216,0) qregless[4];
cx qregless[0],qregless[4];
cx qregless[1],qregless[4];
cx qregless[2],qregless[4];
cx qregless[3],qregless[4];
cx qregless[5],qregless[9];
u3(1.541594019662195,2.7913723796959733,-pi) qregless[5];
cx qregless[0],qregless[5];
cx qregless[1],qregless[5];
cx qregless[2],qregless[5];
cx qregless[3],qregless[5];
cx qregless[4],qregless[5];
cx qregless[6],qregless[9];
u3(1.8595822481541888,-0.0777514342798824,0) qregless[6];
cx qregless[0],qregless[6];
cx qregless[1],qregless[6];
cx qregless[2],qregless[6];
cx qregless[3],qregless[6];
cx qregless[4],qregless[6];
cx qregless[5],qregless[6];
cx qregless[7],qregless[9];
u3(0.7292486063206927,-0.2745466276846109,-pi) qregless[7];
cx qregless[0],qregless[7];
cx qregless[1],qregless[7];
cx qregless[2],qregless[7];
cx qregless[3],qregless[7];
cx qregless[4],qregless[7];
cx qregless[5],qregless[7];
cx qregless[6],qregless[7];
cx qregless[8],qregless[9];
u3(2.0452499401435484,-1.0925023928214674,0) qregless[8];
cx qregless[0],qregless[8];
cx qregless[1],qregless[8];
cx qregless[2],qregless[8];
cx qregless[3],qregless[8];
cx qregless[4],qregless[8];
cx qregless[5],qregless[8];
cx qregless[6],qregless[8];
cx qregless[7],qregless[8];
u3(1.0368254640000032,1.579430266654784,0) qregless[9];
cx qregless[0],qregless[9];
u3(2.5297885440896417,-2.890521540662405,-pi) qregless[0];
cx qregless[1],qregless[9];
u3(0.6105260558088237,-0.4430137085195698,-pi) qregless[1];
cx qregless[2],qregless[9];
u3(2.924458657431965,-1.1628487595917854,-pi) qregless[2];
cx qregless[3],qregless[9];
u3(2.5748407749922118,0.8575991446821432,-pi) qregless[3];
cx qregless[4],qregless[9];
u3(0.24681462267060233,2.176163324566387,0) qregless[4];
cx qregless[5],qregless[9];
u3(2.2442391778450834,0.2707886752855062,0) qregless[5];
cx qregless[6],qregless[9];
u3(0.5002237983271179,-0.7545152110842537,0) qregless[6];
cx qregless[7],qregless[9];
u3(1.9192612708638719,-1.4876032641952888,0) qregless[7];
cx qregless[8],qregless[9];
u3(2.0779707218466736,-0.7659413833270339,0) qregless[8];
u3(1.421066159778744,-0.5183053354057825,-pi) qregless[9];
barrier qregless[0],qregless[1],qregless[2],qregless[3],qregless[4],qregless[5],qregless[6],qregless[7],qregless[8],qregless[9];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
measure qregless[4] -> cregless[4];
measure qregless[5] -> cregless[5];
measure qregless[6] -> cregless[6];
measure qregless[7] -> cregless[7];
measure qregless[8] -> cregless[8];
measure qregless[9] -> cregless[9];
