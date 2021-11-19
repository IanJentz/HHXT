function [rid_Hh,rid_Ch,rid_Cv,rid_solid,C_in,C_out,H_in,H_out] = HXGeom_halfH(model,Wc,Lc,Ws,Ls)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xp = 0.5*Lc*[-1,-1,-1,1,1,1]+Ls*[-1,0,0,0,0,1]+0.5*Wc*[0,0,1,-1,0,0];
yp = 0.5*Wc*[0,1,1]+Ws*[0,0,1];

R1 = [3,4,xp(1),xp(6),xp(6),xp(1),yp(1),yp(1),yp(2),yp(2)]';
R2 = [3,4,xp(1),xp(6),xp(6),xp(1),yp(2),yp(2),yp(3),yp(3)]';
R3 = [3,4,xp(1),xp(2),xp(2),xp(1),yp(1),yp(1),yp(2),yp(2)]';
R4 = [3,4,xp(6),xp(5),xp(5),xp(6),yp(1),yp(1),yp(2),yp(2)]';
P1 = [2,4,xp(2),xp(2),xp(3),xp(3),yp(1),yp(3),yp(3),yp(2)]';
P2 = [2,4,xp(5),xp(5),xp(4),xp(4),yp(1),yp(3),yp(3),yp(2)]';

gd = [R1,R2,R3,R4,P1,P2];
ns = char({'R1','R2','R3','R4','P1','P2'})';
sf = 'R1+R2+R3+R4+P1+P2';

[dl,bt] = decsg(gd,sf,ns);

% build the geometry
G = geometryFromEdges(model,dl);

rid_Hh = [5,2,4,10,7];
rid_Ch = 4;
rid_Cv = [1,2,10,9];
rid_solid = [3,6,8];
C_in = 19;
C_out = 21;
H_in = 12;
H_out = 5;


end

