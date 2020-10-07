function [v1, v2] = callValueFuncOnEdge( self, bc, femEdges, f )
%callValueFuncOnEdge Call user-defined function for the value vector
% The vector of values in a BC of type Value can be defined in a user-written
% function.
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.
numFemEdges = length(femEdges);
[nxe, nye] = self.calcEdgeNormals(femEdges);
if(bc.isVectorized)
% If the user has defined a vectorized function, call it once.
  n1 = self.edges(1,femEdges); n2 = self.edges(2,femEdges);
  loc.x = [self.points(1,n1) self.points(1,n2)];
  loc.y = [self.points(2,n1) self.points(2,n2)];
  loc.nx = [nxe nxe]; loc.ny = [nye nye];
  loc.geomID = bc.appRegionID;
  state.u = [self.uN(:,n1); self.uN(:,n2)];
  state.time = self.time;
  if nargin(f) == 3
    rh = f(self.myPde, loc, state);
  else
    rh = f(loc, state);    
  end
  drh = size(rh);
  ne2 = 2*numFemEdges;
  if(drh(end) ~= ne2)
    error(message('pde:pde2DBCImpl:invalidNumDirFunc', func2str(f), ne2, drh(end)));
  end
  reshape(rh, [], ne2);
  v1 = rh(:, 1:numFemEdges);
  v2 = rh(:, numFemEdges+1:ne2);
else
  firstCall = true;
  % If the user's function is not vectorized, call it for each edge.
  for i=1:numFemEdges
    e = femEdges(i);
    for j=1:2
      n1 = self.edges(j,e);
      loc.x = self.points(1,n1);
      loc.y = self.points(2,n1);
      loc.nx = nxe; loc.ny = nye;
      loc.geomID = bc.appRegionID;
      state.u = self.uN(:,n1);
      state.time = self.time;
      if nargin(f) == 3
        rh = f(self.myPde, loc, state); 
      else
        rh = f(loc, state);   
      end
      rh = rh(:);
      if(firstCall)
        numRows = size(rh,1);
        v1 = zeros(numRows, numFemEdges); v2 = zeros(numRows, numFemEdges);
        firstCall = false;
      end
      if(j == 1)
        v1(:, i) = rh;
      else
        v2(:, i) = rh;
      end
    end
  end
end

end

