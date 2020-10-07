function bc = callNeumannFuncOnFace(self, bc, xyz, sPts, func, faceNormVec, ...
  uPoints, t12)
%callNeumannFuncOnFace Call a user-defined function for Neumann BC
% This undocumented function may be changed or removed in a future release.
% This function handles the case when a function handle has been defined
% for the 'g' or 'h' parameter in a pdeBoundaryConditions entity.

%       Copyright 2014 The MathWorks, Inc.

if(t12 == 1)
  numRows = self.N;
else
  numRows = self.N^2;
end
numFaces = size(xyz, 3);
numIntPts = size(sPts, 2);
xyzPts = zeros(3, numIntPts, numFaces);
if(self.gotTime)
  state.time = self.time;
end
for i=1:numFaces
  xyzPts(:,:,i) = xyz(:,:,i)*sPts;
end
xyzPts = reshape(xyzPts, 3, []);
uPoints = reshape(uPoints, self.N, []);
faceNormVec = reshape(faceNormVec, 3, []);
n = size(xyzPts,2);
if(bc.isVectorized && n>1)
% If the user's function can handle vector arguments and
% we need to evaluate at more than one point, call the function
% once.
  appRegion = self.applicationRegion(xyzPts, faceNormVec);
  state.u = uPoints;
  if nargin(func) == 3
    bc = func(self.myPde, appRegion, state);
  else
    bc = func(appRegion, state);  
  end
  if(t12 == 1)
    self.checkFuncEvalVecN('Neumann', func, bc, n);
  else
    self.checkFuncEvalVecNxN('Neumann', func, bc, n);
  end
  bc = reshape(bc, [], n);
else
  bc = zeros(numRows, n);
  for i=1:n
    appRegion = self.applicationRegion(xyzPts(:,i), faceNormVec(:,i));
    state.u = uPoints(:,i);
    if nargin(func) == 3
        bci = func(self.myPde, appRegion, state);
    else
        bci = func(appRegion, state);
    end
    bci = bci(:);
    if(length(bci) ~= numRows)
      error(message('pde:pde2DBCImpl:invalidNumEntriesDirFunc', func2str(func), numRows));
    end
    bc(:,i) = bci;
  end
end

end

