function [rhV1, rhV2] = callDirichletFuncOnEdge( self, bc, femEdges, f, t12)
%callDirichletFuncOnEdge Call user function for r or h Dirichlet BC terms
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.

if(t12 == 1)
  numRows = self.N;
else
  numRows = self.N^2;
end
numFemEdges = length(femEdges);
    
if(bc.isVectorized)
  n1 = self.edges(1,femEdges); n2 = self.edges(2,femEdges);
  % Include x and y values from both ends of the edge
  loc.x = [self.points(1,n1) self.points(1,n2)];
  loc.y = [self.points(2,n1) self.points(2,n2)];
  loc.geomID = bc.appRegionID;
  state.u = [self.uN(:,n1); self.uN(:,n2)];
  state.time = self.time;
  if nargin(f) == 3
    rh = f(self.myPde, loc, state);
  else
    rh = f(loc, state);  
  end
   
  ne2 = 2*numFemEdges;
  if(t12 == 1)
    self.checkFuncEvalVecN('Dirichlet', f, rh, ne2);
  else
    self.checkFuncEvalVecNxN('Dirichlet', f, rh, ne2);
  end
  rh = reshape(rh, [], ne2);
  rhV1 = rh(:, 1:numFemEdges);
  rhV2 = rh(:, numFemEdges+1:ne2);
else
  rhV1 = zeros(numRows, numFemEdges); rhV2 = zeros(numRows, numFemEdges);
  for i=1:numFemEdges
    e = femEdges(i);
    % Request the BC data at the first node on the edge, then the second node
    for j=1:2
      n1 = self.edges(j,e);
      loc.x = self.points(1,n1);
      loc.y = self.points(2,n1);
      loc.geomID = bc.appRegionID;
      state.u = self.uN(:,n1);
      state.time = self.time;      
      if nargin(f) == 3
         rh = f(self.myPde, loc, state);
      else
         rh = f(loc, state);  
      end   
      rh = rh(:);
      if(length(rh) ~= numRows)
        error(message('pde:pde2DBCImpl:invalidNumEntriesDirFunc', func2str(f), numRows));
      end
      if(j == 1)
        rhV1(:, i) = rh;
      else
        rhV2(:, i) = rh;
      end
    end
  end
end

end

