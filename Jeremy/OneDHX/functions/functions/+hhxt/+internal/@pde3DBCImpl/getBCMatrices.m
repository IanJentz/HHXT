function [Q,G,H,R] = getBCMatrices(self,p,meshAssoc,u,time)
  %getBCMatrices Calculate global matrices for 3D boundary conditions
  %
  % This function iterates over all BCs applied to the faces of a 3D model.
  % First it gets the element faces for the particular geometry face.
  % Then it calls the appropriate lower-level function to compute the BC
  % terms for each element face.
  % This undocumented function may be changed or removed in a future release.
  
  %       Copyright 2014 The MathWorks, Inc.

import pde.internal.*;

self.points = p;
self.numPoints = size(p,2);
neqn = self.N*self.numPoints;

self.gotU = nargin > 3 && ~isempty(u);
self.gotTime = nargin > 4 && ~isempty(time);
if(self.gotU)
  self.uN = reshape(u,[],self.N)';
else
  self.uN = zeros(self.N, self.numPoints);
end
if(self.gotTime)
  self.time = time;
else
  self.time = [];
end


Q = sparse(neqn,neqn);
G = sparse(neqn,1);

% map of Dirichlet BCs with node ID as the key
dirBoundaryConditions = containers.Map('KeyType','uint32','ValueType','any');

% map of neumann BCs with geom face ID as key. 
% This "map" just tracks which faces already have Neumann BCs applied to
% them so they are applied only once.
neuBoundaryConditions = containers.Map('KeyType','uint32','ValueType','logical');

% iterate, last to first so only last Neumann BC on a face is applied
for i=length(self.faceBoundaryConditions):-1:1 
  bcsi = self.faceBoundaryConditions(i);
  bcsi.elemFaces = meshAssoc.getElementFaces(bcsi.appRegionID);
  % Call the appropriate BC function for each type
  if(bcsi.bcType == pdeBCImpl.valueBC)
    bcsi.nodes = self.getNodesForFaces(bcsi.elemFaces);
    bcsi.r = bcsi.term1;
    bcsi.h = bcsi.term2;
    self.setValueBCOnFace(bcsi, dirBoundaryConditions);
  elseif(bcsi.bcType == pdeBCImpl.neumannBC)
    % Last Neumann BC on a face is the right one. So, insert only the last one.
    if ~neuBoundaryConditions.isKey(bcsi.faceID)
      neuBoundaryConditions(bcsi.faceID) = true;
      bcsi.g = bcsi.term1;
      bcsi.q = bcsi.term2;
      [Qi, Gi] = setNeumannBCOnFace(self, bcsi);
      Q = Q + Qi;
      G = G + Gi;
    end
  elseif(bcsi.bcType == pdeBCImpl.dirichletBC)
    bcsi.nodes = self.getNodesForFaces(bcsi.elemFaces);
    bcsi.r = bcsi.term1;
    bcsi.h = bcsi.term2;
    setDirichletBCOnFace(self, bcsi, dirBoundaryConditions);
  else
    error(message('pde:pde2DBCImpl:invalidBCType', bcsi.bcType));
  end
end

%
% Create the H and R matrices
%
% preallocate some space to hold the triplets
numDirBCs = length(dirBoundaryConditions);
preAllocSize = max(numDirBCs*self.N, 1);
tripR = zeros(2, preAllocSize);
tripH = zeros(3, preAllocSize);
bcNodes = dirBoundaryConditions.keys;
nodalVals = dirBoundaryConditions.values;
% actual counts of constraints and terms
consCount = 0; hCount = 0; rCount = 0;
for i=1:length(dirBoundaryConditions)
  nodeID = bcNodes{i};
  hr = nodalVals{i};
  hNodeI = hr.h;
  rNodeI = hr.r;
  for r=1:self.N
    if(any(hNodeI(r,:)~=0))
      rr = rNodeI(r);
      consCount = consCount + 1;
      if(rr~=0)
        rCount = rCount + 1;
        tripR(1, rCount) = consCount;
        tripR(2, rCount) = rr;
      end
      for c=1:self.N
        hrc = hNodeI(r,c);
        if(hrc~=0)
          eq = (c-1)*self.numPoints + nodeID;
          hCount = hCount + 1;
          tripH(1, hCount) = consCount;
          tripH(2, hCount) = eq;
          tripH(3, hCount) = hrc;
          hNodeI(r,c) = hrc;
        end
      end
    end
  end
end


if(consCount)
  H = sparse(tripH(1, 1:hCount), tripH(2, 1:hCount), tripH(3, 1:hCount), consCount, neqn);
  R = sparse(tripR(1, 1:rCount), 1, tripR(2, 1:rCount), consCount, 1);
else
  H = sparse(1,neqn);
  R = sparse(neqn,1);
end

end


