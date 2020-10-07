function [G] = Geom_Rectangle(model,L,W)
% Geom_Rectangle - Build a rectangular pde geometry using the
% geometryFromEdges() function.
%   G = Geom_Rectangle(model,L,W) contructs a 2D rectangular geometry of width
%   W and length L, both in meters.

% See also pde.geometryFromEdges, DECSG, PDEGEOM, pde.AnalyticGeometry

    % a sketch of the rectangular geometry
    %  
    % (x1,y2)           (x2,y2)
    %    +-----------------+
    %    |                 |
    %    |        .(0,0)   | W
    %    |                 |
    %    +-----------------+
    % (x1,y1)     L     (x2,y1)
    
    % define the x and y locations of the rectangles corners
    x1 = -0.5*L;
    x2 =  0.5*L;
    y1 = -0.5*W;
    y2 =  0.5*W;
    
    % define a Rectangle 
    R1 = [3,4,x1,x2,x2,x1,y2,y2,y1,y1]';
    
    % assemble inputs needed for decsg() (see decsg in help browser for details
    % on input syntax)
    gd = R1;
    ns = char('R1')';
    sf = 'R1';

    % Decompose constructive solid 2-D geometry into minimal regions
    [dl,bt] = decsg(gd,sf,ns);
    
    % create a geometry from the edges generated in decsg
    G = geometryFromEdges(model,dl);

end

