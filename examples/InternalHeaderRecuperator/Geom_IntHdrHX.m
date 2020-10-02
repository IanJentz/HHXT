function G = IntHdrHX(model)
%IntHdrHX - Builds the internally headered heat exchanger geometry using
%the geometryFromEdges() function.
% See also pde.geometryFromEdges, DECSG, PDEGEOM, pde.AnalyticGeometry

% use decsg to build the geometry out intersecting rectangles and
% circles 
x1 = -0.5*30.125*25.4e-3;
x2 = x1 + 1.1964*25.4e-3;
x3 = x1 + 3.0723*25.4e-3;
x4 = 0;
x7 = -x1;
x6 = -x2;
x5 = -x3;

y1 = -0.5*4.2580*25.4e-3;
y2 = y1 + 0.4403*25.4e-3;
y3 = y1 + 1.1702*25.4e-3;
y4 = 0;
y7 = -y1;
y6 = -y2;
y5 = -y3;

r1 = 0.5512*25.4e-3;

x01 = x1 + 1.0609*25.4e-3;
x02 = x1 + 1.5683*25.4e-3;
x03 = x1 + 1.7242*25.4e-3;
x04 = x1 + 2.8131*25.4e-3;
x05 = x1 + 2.9132*25.4e-3;

y01 = y1 + 1.4362*25.4e-3;

zbuff = zeros(1,2);
R1 = [3,4,x1,-x1,-x1,x1,y1,y1,-y1,-y1,zbuff]';
R2 = [3,4,x3,-x3,-x3,x3,y2,y2,-y2,-y2,zbuff]';
R3 = [3,4,x05,x3,x3,x05,y2,y2,-y2,-y2,zbuff]';
R4 = [3,4,-x05,-x3,-x3,-x05,y2,y2,-y2,-y2,zbuff]';

zbuff = zeros(1,8);
C1 = [1,x2,y3,r1,zbuff]';
C2 = [1,-x2,y3,r1,zbuff]';
C3 = [1,-x2,-y3,r1,zbuff]';
C4 = [1,x2,-y3,r1,zbuff]';

zbuff = zeros(1,2);
P1 = [2,4,x01,x02,x04,x03,y01,y2,y2,y4,zbuff]';
P2 = [2,4,-x01,-x02,-x04,-x03,y01,y2,y2,y4,zbuff]';
P3 = [2,4,-x01,-x02,-x04,-x03,-y01,-y2,-y2,-y4,zbuff]';
P4 = [2,4,x01,x02,x04,x03,-y01,-y2,-y2,-y4,zbuff]';
P5 = [2,5,x03,x04,x3,x3,x04,y4,y2,y2,-y2,-y2]';
P6 = [2,5,-x03,-x04,-x05,-x05,-x04,y4,y2,y2,-y2,-y2]';

gd = [R1,R2,R3,R4,C1,C2,C3,C4,P1,P2,P3,P4,P5,P6];
ns = char('R1','R2','R3','R4','C1','C2','C3','C4','P1','P2','P3','P4','P5','P6')';
sf = ['(R1-C1-C2-C3-C4)+R2+R3+R4+P5+P6',...
        '+(P1-C1)+(P2-C2)+(P3-C3)+(P4-C4)'];

[dl,bt] = decsg(gd,sf,ns);
  
G = geometryFromEdges(model,dl);

end

