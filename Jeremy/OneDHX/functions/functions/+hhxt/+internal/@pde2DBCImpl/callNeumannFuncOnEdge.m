function qg = callNeumannFuncOnEdge( self, bc, femEdges, f, t12 )
%callNeumannFuncOnEdge Call user function for g or q Neumann BC terms
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.

if(t12 == 1)
  numRows = self.N;
else
  numRows = self.N^2;
end
numFemEdges = length(femEdges);
[nxe, nye] = self.calcEdgeNormals(femEdges);
if(bc.isVectorized)
  % If the user has defined a vectorized function, call it once.
  n1 = self.edges(1,femEdges); n2 = self.edges(2,femEdges);
  loc.x = (self.points(1,n1) + self.points(1,n2))/2;
  loc.y = (self.points(2,n1) + self.points(2,n2))/2;
  loc.nx = nxe; loc.ny = nye;
  state.u = (self.uN(:,n1) + self.uN(:,n2))/2;
  state.time = self.time;
  loc.geomID = bc.appRegionID;
  if nargin(f) == 3
    qg = f(self.myPde, loc, state);
  else
    qg = f(loc, state);  
  end  
  if(t12 == 1)
    self.checkFuncEvalVecN('Neumann', f, qg, numFemEdges);
  else
    self.checkFuncEvalVecNxN('Neumann', f, qg, numFemEdges);
  end
else
  qg = zeros(numRows, numFemEdges);
  % If the user's function is not vectorized, call it for each edge.
  for i=1:numFemEdges
    e = femEdges(i);
    n1 = self.edges(1,e); n2 = self.edges(2,e);
    loc.x = (self.points(1,n1) + self.points(1,n2))/2;
    loc.y = (self.points(2,n1) + self.points(2,n2))/2;
    loc.nx = nxe(i); loc.ny = nye(i);
    state.u = (self.uN(:,n1) + self.uN(:,n2))/2;
    state.time = self.time;
    loc.geomID = bc.appRegionID;
    if nargin(f) == 3
        qgi = f(self.myPde, loc, state); 
    else
        qgi = f(loc, state); 
    end
    qgi = qgi(:);
    if(length(qgi) ~= numRows)
      error(message('pde:pde2DBCImpl:invalidNumEntriesDirFunc', func2str(f), numRows));
    end
    qg(:, i) = qgi;
  end
end

end

