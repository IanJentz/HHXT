function [nxe, nye] = calcEdgeNormals( self, femEdges )
%calcEdgeNormals Calculate the unit normals to a collection of 2D element edges
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.
e = self.edges(:,femEdges);
ne = size(e, 2);
% Coordinates of all points on femEdges
xp = [self.points(1,e(1,:)) self.points(1,e(2,:))];
yp = [self.points(2,e(1,:)) self.points(2,e(2,:))];
ls = find(e(6,:)==0 & e(7,:)>0); % External region to the left
rs = find(e(7,:)==0 & e(6,:)>0); % External region to the right
% any internal edges will have normals equal NaN
nxe = NaN*zeros(1,ne);
nye = NaN*zeros(1,ne);
dx = xp(ne+1:2*ne)-xp(1:ne);
dy = yp(ne+1:2*ne)-yp(1:ne);
dl = sqrt(dx.^2+dy.^2);
% normals for edges with external region to the right
nxe(rs) = dy(rs)./dl(rs);
nye(rs) = -dx(rs)./dl(rs);
% normals for edges with external region to the left
nxe(ls) = -dy(ls)./dl(ls);
nye(ls) = dx(ls)./dl(ls);
end

