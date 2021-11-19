function [rid_Hh,rid_Ch,rid_Cv,rid_solid,C_in,C_out,H_in,H_out] = HXGeom_S(model,Wc,Lc,Ws,Ls)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xp = 0.5*Lc*[-1,-1,-1,1,1,1]+Ls*[-1,0,0,0,0,1]+Wc*[0,0,1,-1,0,0];
yp = 0.5*Wc*[-1,-1,1,1]+Ws*[-1,0,0,1];

R1 = [3,4,xp(1),xp(6),xp(6),xp(1),yp(1),yp(1),yp(4),yp(4)]';
R2 = [3,4,xp(1),xp(6),xp(6),xp(1),yp(1),yp(1),yp(2),yp(2)]';
R3 = [3,4,xp(1),xp(6),xp(6),xp(1),yp(3),yp(3),yp(4),yp(4)]';
R4 = [3,4,xp(1),xp(1),xp(2),xp(2),yp(2),yp(3),yp(3),yp(2)]';
R5 = [3,4,xp(6),xp(6),xp(5),xp(5),yp(2),yp(3),yp(3),yp(2)]';
P1 = [2,4,xp(2),xp(2),xp(3),xp(3),yp(4),yp(2),yp(3),yp(4)]';
P2 = [2,4,xp(1),xp(1),xp(2),xp(3),yp(3),yp(2),yp(2),yp(3)]';
P3 = [2,4,xp(5),xp(5),xp(4),xp(4),yp(1),yp(3),yp(2),yp(1)]';
P4 = [2,4,xp(6),xp(6),xp(5),xp(4),yp(2),yp(3),yp(3),yp(2)]';


gd = [R1,R2,R3,R4,R5,P1,P2,P3,P4];
ns = char({'R1','R2','R3','R4','R5','P1','P2','P3','P4'})';
sf = 'R1+R2+R3+R4+R5+P1+P2+P3+P4';

[dl,bt] = decsg(gd,sf,ns);

% build the geometry
G = geometryFromEdges(model,dl);

rid_Hh = [6,3,7,10,9];
rid_Ch = 7;
rid_Cv = [4,3,10,11];
rid_solid = [5,1,8,2];
C_in = 6;
C_out = 15;
H_in = 18;
H_out = 9;


end

