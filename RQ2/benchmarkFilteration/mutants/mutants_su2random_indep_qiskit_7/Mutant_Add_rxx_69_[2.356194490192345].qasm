OPENQASM 2.0;
include "qelib1.inc";
qreg qregless[7];
creg cregless[7];
u3(1.4368347742816614,1.6369627429575306,-pi) qregless[0];
u3(0.13038834331032634,1.0625547235745714,0) qregless[1];
cx qregless[0],qregless[1];
u3(2.301856027570534,-2.5865372311326773,-pi) qregless[2];
cx qregless[0],qregless[2];
cx qregless[1],qregless[2];
u3(1.5783117544539513,1.1646500873100205,-pi) qregless[3];
cx qregless[0],qregless[3];
cx qregless[1],qregless[3];
cx qregless[2],qregless[3];
u3(3.1322119352256292,-0.2928382424047804,0) qregless[4];
cx qregless[0],qregless[4];
cx qregless[1],qregless[4];
cx qregless[2],qregless[4];
cx qregless[3],qregless[4];
u3(1.4124389803026796,0.02480768898038388,0) qregless[5];
cx qregless[0],qregless[5];
cx qregless[1],qregless[5];
cx qregless[2],qregless[5];
cx qregless[3],qregless[5];
cx qregless[4],qregless[5];
u3(1.2444656817555668,-3.064986403423018,0) qregless[6];
cx qregless[0],qregless[6];
u3(1.1773372206208805,-2.248311899378857,-pi) qregless[0];
cx qregless[1],qregless[6];
u3(2.434570523812673,-0.7958234754631421,-pi) qregless[1];
cx qregless[0],qregless[1];
cx qregless[2],qregless[6];
u3(1.748262901313365,1.0941137716709264,-pi) qregless[2];
cx qregless[0],qregless[2];
cx qregless[1],qregless[2];
cx qregless[3],qregless[6];
u3(1.8339114230470699,2.776119709759085,0) qregless[3];
cx qregless[0],qregless[3];
cx qregless[1],qregless[3];
cx qregless[2],qregless[3];
cx qregless[4],qregless[6];
u3(0.5166404252966228,-0.4146023075677032,-pi) qregless[4];
cx qregless[0],qregless[4];
cx qregless[1],qregless[4];
cx qregless[2],qregless[4];
cx qregless[3],qregless[4];
cx qregless[5],qregless[6];
u3(1.7933732440688746,0.7399517487893483,-pi) qregless[5];
cx qregless[0],qregless[5];
cx qregless[1],qregless[5];
cx qregless[2],qregless[5];
cx qregless[3],qregless[5];
cx qregless[4],qregless[5];
u3(2.8742785055981956,0.0825500125799099,-pi) qregless[6];
cx qregless[0],qregless[6];
u3(2.1966192898367836,-1.2522384758651306,-pi) qregless[0];
cx qregless[1],qregless[6];
u3(2.5067461861055573,-2.4254077858804965,-pi) qregless[1];
cx qregless[0],qregless[1];
cx qregless[2],qregless[6];
u3(1.223818747839898,2.0651656802006935,-pi) qregless[2];
cx qregless[0],qregless[2];
cx qregless[1],qregless[2];
cx qregless[3],qregless[6];
u3(3.005579583727834,-2.846934388642458,-pi) qregless[3];
cx qregless[0],qregless[3];
cx qregless[1],qregless[3];
cx qregless[2],qregless[3];
cx qregless[4],qregless[6];
u3(0.5739760098973872,0.7934855547557511,-pi) qregless[4];
cx qregless[0],qregless[4];
cx qregless[1],qregless[4];
cx qregless[2],qregless[4];
rxx(3*pi/4) qregless[3],qregless[4];
cx qregless[3],qregless[4];
cx qregless[5],qregless[6];
u3(2.0058195038543025,-2.8426000178928454,0) qregless[5];
cx qregless[0],qregless[5];
cx qregless[1],qregless[5];
cx qregless[2],qregless[5];
cx qregless[3],qregless[5];
cx qregless[4],qregless[5];
u3(0.5683728542359917,-1.1354532936221062,0) qregless[6];
cx qregless[0],qregless[6];
u3(1.2500242582094414,1.0368254640000032,0) qregless[0];
cx qregless[1],qregless[6];
u3(0.8994360763247731,-0.6752586753862846,-pi) qregless[1];
cx qregless[2],qregless[6];
u3(2.2094986973106154,0.5872288522304125,0) qregless[2];
cx qregless[3],qregless[6];
u3(1.541594019662195,2.0175663513732243,-pi) qregless[3];
cx qregless[4],qregless[6];
u3(1.8595822481541888,0.9497161489686778,0) qregless[4];
cx qregless[5],qregless[6];
u3(0.7292486063206927,-0.7281303932915781,-pi) qregless[5];
u3(2.045249940143549,-0.3502202738938198,0) qregless[6];
barrier qregless[0],qregless[1],qregless[2],qregless[3],qregless[4],qregless[5],qregless[6];
measure qregless[0] -> cregless[0];
measure qregless[1] -> cregless[1];
measure qregless[2] -> cregless[2];
measure qregless[3] -> cregless[3];
measure qregless[4] -> cregless[4];
measure qregless[5] -> cregless[5];
measure qregless[6] -> cregless[6];
