function rMat = callValueFuncOnFace(self, bc, nodes, func, faceNormals)
  %callValueFuncOnFace Call a user-defined function for value type BC
  %
  % This undocumented function may be changed or removed in a future release.
  % This function handles the case when a function handle has been defined
  % for the 'u' parameter in a pdeBoundaryConditions entity.
  
  %       Copyright 2014 The MathWorks, Inc.

n = length(nodes);
p = self.points;
if(self.gotTime)
  state.time = self.time;
end
if(bc.isVectorized)
  appRegion = self.applicationRegion(p(:,nodes), faceNormals);
  state.u = self.uN(:,nodes);
  if nargin(func) == 3
    rMat = func(self.myPde, appRegion, state);
  else
    rMat = func(appRegion, state);  
  end
  numRet = size(rMat, ndims(rMat));
  if(numRet ~= n)
    error(message('pde:pde2DBCImpl:invalidNumDirFunc', func2str(func), n, numRet));
  end
else
  % If the user's function is not vectorized, call it for each face.
  xyzPts = p(:,nodes);
  for i=1:n
    ni = nodes(i);
    appRegion = self.applicationRegion(xyzPts(:,i), faceNormals(:,i));
    state.u = self.uN(:,ni);
    if nargin(func) == 3
        bci = func(self.myPde, appRegion, state);
    else
        bci = func(appRegion, state);
    end
    if(i == 1)
      rMat = zeros(length(bci), n);
    end
    rMat(:,i) = bci(:);
  end
end

end

