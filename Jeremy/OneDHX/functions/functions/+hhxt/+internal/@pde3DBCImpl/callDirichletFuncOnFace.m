function rh = callDirichletFuncOnFace( self, bc, nodes, func, faceNormals, t12 )
%callDirichletFuncOnFace Call a user-defined function for a Dirichlet BC
% This undocumented function may be changed or removed in a future release.
% This function handles the case when a function handle has been defined
% for the 'h' or 'r' parameter in a pdeBoundaryConditions entity.

%       Copyright 2014-2015 The MathWorks, Inc.

% If t12 is 1, this function is being called for an r-vector, otherwise
% and h-matrix
if(t12 == 1)
  numRows = self.N;
else
  numRows = self.N^2;
end
n = length(nodes);
p = self.points;
if(self.gotTime)
  state.time = self.time;
end
if(bc.isVectorized)
  appRegion = self.applicationRegion(p(:,nodes), faceNormals);
  state.u = self.uN(:,nodes);
  if nargin(func) == 3
    rh = func(self.myPde, appRegion, state);
  else
    rh = func(appRegion, state);  
  end
  if(t12 == 1)
    self.checkFuncEvalVecN('Dirichlet', func, rh, n);
  else
    self.checkFuncEvalVecNxN('Dirichlet', func, rh, n);
  end
  rh = reshape(rh, [], n);
else
  xyzPts = p(:,nodes);
  rh = zeros(numRows, n);
  for i=1:n
    ni = nodes(i);
    appRegion = self.applicationRegion(xyzPts(:,i), faceNormals(:,i));
    state.u = self.uN(:,ni);
    if nargin(func) == 3
        bci = func(self.myPde, appRegion, state); 
        bci = bci(:);
    else
        bci = func(appRegion, state); 
        bci = bci(:);
    end
    if(length(bci) ~= numRows)
      error(message('pde:pde2DBCImpl:invalidNumEntriesDirFunc', func2str(func), numRows));
    end
    rh(:,i) = bci(:);
  end
end

end

