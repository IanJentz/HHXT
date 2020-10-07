function [ Qi, Gi ] = setNeumannBCOnFace( self, bci )
%setNeumannBCOnFace Calculate g and q entries for a Neumann BC
% This undocumented function may be changed or removed in a future release.
% This function handles the case when a 'g' and/or 'q' parameter
% has been included in a pdeBoundaryConditions entity.

%       Copyright 2014 The MathWorks, Inc.


import pde.internal.elements.*;

N2 = self.N*self.N;
neqn = self.N*self.numPoints;
faceID = bci.appRegionID;

isFuncQ = isa(bci.q, 'function_handle');
isFuncG = isa(bci.g, 'function_handle');
anyFunc = isFuncQ || isFuncG;
anyQ = isFuncQ || any(bci.q(:)); anyG = isFuncG || any(bci.g(:));

if(~ (anyQ || anyG))
  Qi = sparse(neqn,neqn);
  Gi = sparse(neqn,1);
  return;
end
[numFaceNodes, numFaces] = size(bci.elemFaces);
% preallocate some space for sparse triplets
preAllocSize = 3*numFaces*self.N;
if(anyQ)
  qTriplets = zeros(3, preAllocSize);
  numQTriplets = 0;
end
if(anyG)
  gTriplets = zeros(2, preAllocSize);
  numGTriplets = 0;
end
if(numFaceNodes == 3)
  sfCalc = @shapeTri3;
  dsfCalc = @dShapeTri3;
  numPoints = 1;
elseif(numFaceNodes == 6)
  sfCalc = @shapeTri6;
  dsfCalc = @dShapeTri6;
  numPoints = 3;
else
  error(message('pde:pde3DBCImpl:illegalFaceDefn', numFaceNodes));
end

% compute shape functions and derivatives at gauss points
intRule=GaussianIntegrationRule(GaussianIntegrationRule.Triangle, numPoints);
intWts = intRule.wts;
intPoints = intRule.points;
dsPts = zeros(numFaceNodes, 2, numPoints);
for n=1:numPoints
  dsPts(:,:,n) = dsfCalc(intPoints(:,n));
end
sPts = sfCalc(intPoints);

if(anyG)
  Ge = zeros(numFaceNodes, self.N);
end

if(anyQ)
  Qe = zeros(numFaceNodes, numFaceNodes, N2);
end

if(~isempty(bci.q) && ~isFuncQ)
    qi = bci.q(:);
    if(any(qi) && length(qi) ~= N2)
        error(message('pde:pde3DBCImpl:invalidLenQ', faceID, N2));
    elseif(~any(qi) && length(qi) ~= N2 && length(qi) ~= 1)
        error(message('pde:pde3DBCImpl:invalidLenQ', faceID, N2));
    end
end
if(~isempty(bci.g) && ~isFuncG)
    gi = bci.g(:);
    if(any(gi) && length(gi) ~= self.N)
        error(message('pde:pde3DBCImpl:invalidLenG', numFaceNodes, self.N));
    elseif(~any(gi) && length(gi) ~= self.N && length(gi) ~= 1)
        error(message('pde:pde3DBCImpl:invalidLenG', numFaceNodes, self.N));
    end
end

%  calculate face normals
xyzAllFaceNodes = self.points(:, bci.elemFaces(:));
xyzAllFaceNodes = reshape(xyzAllFaceNodes, 3, numFaceNodes, numFaces);
faceNorm = zeros(numPoints, numFaces);
if(anyFunc)
  uPoints = zeros(self.N, numPoints, numFaces);
  faceNormVec = zeros(3, numPoints, numFaces);
end
for i=1:numFaces
  xyz = xyzAllFaceNodes(:,:,i);
  if(anyFunc)
    faceNodes = bci.elemFaces(:,i);
    uFaceNodes = self.uN(:, faceNodes);
    uPoints(:,:,i) = uFaceNodes*sPts;
  end
  for j=1:numPoints
    dxyzDrs = xyz*dsPts(:,:,j);
    v3 = cross3(dxyzDrs(:,1),dxyzDrs(:,2));
    faceNorm(j,i) = norm(v3);
    if(anyFunc)
      faceNormVec(:,j,i) = v3/faceNorm(j,i);
    end
  end
end
if(isFuncQ)
  faceQ = self.callNeumannFuncOnFace(bci,xyzAllFaceNodes, sPts, bci.q, ...,
    faceNormVec, uPoints, 2);
end
if(isFuncG)
  faceG = self.callNeumannFuncOnFace(bci,xyzAllFaceNodes, sPts, bci.g, ...
    faceNormVec, uPoints, 1);
end

fCol = 1;
% Iterate over all element faces in this BC
for i=1:numFaces
  coni = bci.elemFaces(:,i);
  % Iterate over the integration points
  for n=1:numPoints 
    shn = sPts(:,n);
    detJ = faceNorm(n, i); 
    wt = intWts(n);
    detWt = detJ*wt;
    if(anyG)
      %bci.g
      if(isFuncG)
        gi = faceG(:,fCol);
      end
      detShp = detWt*shn;
      for j=1:self.N
        Ge(:,j) = Ge(:,j) + gi(j)*detShp;
      end
    end
    if(anyQ)
      if(isFuncQ)
        qi = faceQ(:,fCol);
      end
      detShpShp = detWt*(shn*shn');
      for j=1:N2
        Qe(:,:,j) = Qe(:,:,j) +qi(j)*detShpShp;
      end
    end
    fCol = fCol + 1;
  end
  if(anyG)
    % copy the face g-vector to triplets
    for j=1:self.N
      offset = (j-1)*self.numPoints;
      gj = Ge(:,j);
      nzInd = gj ~= 0;
      nnz = sum(nzInd);
      gTriplets(:, numGTriplets+1:numGTriplets+nnz) = [coni(nzInd)'+offset; gj(nzInd)'];
      numGTriplets = numGTriplets + nnz;
    end
    % set Ge back to zero for the next face
    Ge(:) = 0;
  end
  if(anyQ)
    % copy the face q-matrix to triplets
    jj = 1;
    for j1=1:self.N
      offset1 = (j1-1)*self.numPoints;
      for j2=1:self.N
        offset2 = (j2-1)*self.numPoints;
        qjk = reshape(Qe(:, :, jj)', [], 1);
        nzInd = qjk ~= 0;
        nnz = sum(nzInd);
        if(nnz)
          gRows = repmat(coni' + offset2, numFaceNodes, 1); gRows = gRows(:);
          gCols = repmat(coni' + offset1, 1, numFaceNodes);
          qTriplets(:,numQTriplets+1:numQTriplets + nnz) = [gRows(nzInd)'; gCols(nzInd); qjk(nzInd)'];
          numQTriplets = numQTriplets + nnz;
        end
        jj = jj + 1;
      end
    end
    % set Qe back to zero for the next face
    Qe(:) = 0;
  end
end

if(anyG && numGTriplets)
  Gi = sparse(gTriplets(1, 1:numGTriplets), 1, gTriplets(2, 1:numGTriplets), ...
    neqn,1);
else
  Gi = sparse(neqn,1);
end
if(anyQ && numQTriplets)
  Qi = sparse(qTriplets(1, 1:numQTriplets), qTriplets(2, 1:numQTriplets), ...
    qTriplets(3, 1:numQTriplets), neqn, neqn);
else
  Qi = sparse(neqn,neqn);
end

end

% This highly specialized version of the cross function with no error
% checking is used to improve performance.
function c=cross3(a,b) 
c=[a(2)*b(3)-a(3)*b(2); 
   a(3)*b(1)-a(1)*b(3); 
   a(1)*b(2)-a(2)*b(1)];
end


