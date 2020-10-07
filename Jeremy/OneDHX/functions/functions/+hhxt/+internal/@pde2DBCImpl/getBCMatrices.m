function [q,g,h,r] = getBCMatrices(self,p,e,u,time)
%getBCMatrices boundary function compatible with pdebound "boundary file" definition
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.

import pde.internal.*;

ne = size(e,2); % number of boundary edges
N = self.N;
% Set all the returned matrices to their default values-- zero-Neumann
% conditions on all edges.
q = zeros(N^2, ne);
g = zeros(N, ne);
h = zeros(N^2, 2*ne);
r = zeros(N,2*ne);
eyeN = eye(N); 
self.gotU = nargin > 3 && ~isempty(u);
self.gotTime = nargin > 4 && ~isempty(time);

self.points = p;
self.numPoints = size(p, 2);
self.edges = e;
self.numEdges = size(e, 2);
self.hDefault = eyeN(:);
if(self.gotU)
  self.uN = reshape(u,[],N)';
else
  self.uN = zeros(N, size(p,2));
end
if(self.gotTime)
  self.time = time;
else
  self.time = [];
end
% allEdgeBCs is a map with the key being the id of the edge
% the boundary condition is associated with
allEdgeBCs = self.edgeBoundaryConditions.keys();
for i=1:length(self.edgeBoundaryConditions)
  gi = allEdgeBCs{i};
  bcsi = self.edgeBoundaryConditions(gi);
  if(~isempty(bcsi))
    femEdgesI = find(e(5,:)==gi);
	  % In general, there may be more than one boundary condition definition on an edge
    % Iterate over all BCs on this edge
    for j = 1:size(bcsi,2)
      bcij = bcsi(j);
      % Call the appropriate BC function for each type
      if(bcij.bcType == pdeBCImpl.valueBC)
        [h, r] = setValueBCOnEdge(self, bcij, femEdgesI, h, r);
      elseif(bcij.bcType == pdeBCImpl.neumannBC)
        [q, g] = setNeumannBCOnEdge(self, bcij, femEdgesI, q, g);
      elseif(bcij.bcType == pdeBCImpl.dirichletBC)
        [h, r] = setDirichletBCOnEdge(self, bcij, femEdgesI, h, r);
      else
        error(message('pde:pde2DBCImpl:invalidBCType', bcij.bcType));
      end
    end
  end
end
end

