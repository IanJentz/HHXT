function G = Geom_SideHeaders(model,L,W,W_head,alpha)
% GEOM_SIDEHEADERS   Build a HHXT geometry with the side header
% configurations shown below:
%   
%   6----------------5---4
%   |\                \  |
%  W| \        . (0,0) \ |
%   |  \                \|
%   1---2----------------3
%   W_head    L


% postions of the 6 points that describe the geometry ^see diagram above
if alpha == 1
    points = [ -0.5*L , -0.5*W ;...
           0.5*L , -0.5*W ;...
           0.5*L , 0.5*W ;...
           -0.5*L , 0.5*W] ;
else 
points = [ -0.5*L , -0.5*W ;...
           -0.5*L+W_head , -0.5*W ;...
           0.5*L , -0.5*W ;...
           0.5*L , 0.5*W ;...
           0.5*L-W_head , 0.5*W ;...
           -0.5*L , 0.5*W] ;
end
x = points(:,1); y = points(:,2);

if alpha == 1
    % main rectangle describing the WxL envelope
    R1 = [3;4;x(1);x(2);x(3);x(4);y(1);y(2);y(3);y(4)];
else
    R1 = [3;4;x(1);x(3);x(4);x(6);y(1);y(3);y(4);y(6)];

    % two triangles that describe the header regions
    % are described using polygons
    P1 = [2;3;x(1);x(2);x(6);y(1);y(2);y(6)];
    P2 = [2;3;x(3);x(4);x(5);y(3);y(4);y(5)];
    
    % note we need to pad the P columns so that they are the same size as R1
    P1 = padarray(P1,length(R1)-length(P1),0,'post'); % pads with 0 below existing P1
    P2 = padarray(P2,length(R1)-length(P2),0,'post'); % pads with 0 below existing P2
end

% assemble inputs needed for decsg() (see decsg in help browser for details
% on input syntax)
if alpha == 1
    gd = [R1];    % matrix of geom features being assembled
    ns = char('R1')'; % variable name of each column in the ^ gd matrix
    sf = 'R1';
else
    gd = [R1,P1,P2];    % matrix of geom features being assembled
    ns = char('R1','P1','P2')'; % variable name of each column in the ^ gd matrix
    sf = 'R1+P1+P2'; % operations used to create the geometry from the geom feature
    % i.e. build a geometry out of the combination of R1 P1 and P2
    % other operations are:
    %  + combination of features, e.g R1+P1 gives the combined area of R1 and P1
    %  * union of features, e.g. R1*P1 gives the overlap of R1 and P1
    %  - subtract feature, e.g. R1-P1 gives area of R1 excluding that overlaped by P1
    %  () order of operations, e.g. (R1-P1)+P2 first subtracts then combines
end
% Decompose constructive solid 2-D geometry into minimal regions
[dl,bt] = decsg(gd,sf,ns);

% create a geometry from the edges generated in decsg
G = geometryFromEdges(model,dl);

end

