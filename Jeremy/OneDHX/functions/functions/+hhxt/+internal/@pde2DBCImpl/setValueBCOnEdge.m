function [h, r] = setValueBCOnEdge(self, bc, femEdges, h, r)
%setValueBCOnEdge Implement type 'Value" BC for a collection of edges
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.

edgeID = bc.appRegionID;
numFemEdges = length(femEdges);
ne = self.numEdges;
% Both the first and second terms are optional. But to get to this point,
% at least one must have been included to determine the type.
%
% Check to see if a second term has been included.
if(~isempty(bc.term2))
  % we have a second term
  if ~isnumeric(bc.term2)
    error(message('pde:pde2DBCImpl:term2NotNum', edgeID));
  end
  t2 = bc.term2;
  rr1 = t2(:);
  lenRr1 = length(rr1);
  if(lenRr1 > self.N)
    eId = 'pde:pde2DBCImpl:t2TooLarge';
    err = message(eId, lenRr1, self.N, edgeID);
    throwAsCaller(MException(eId, err.getString()));
  end
  hN = getHInd(self.N, t2, bc.appRegionID);
  hAllEdges = ones(lenRr1, numFemEdges);
  h(hN,femEdges) = hAllEdges;
  h(hN,femEdges+ne) = hAllEdges;
else
  % no second term, default h-matrix (all dofs)
  rr1 = 1:self.N;
  lenRr1 = self.N;
  hAllEdges = repmat(self.hDefault, 1, numFemEdges);
  h(:,femEdges) = hAllEdges;
  h(:,femEdges+ne) = hAllEdges;
end
%
% Now handle the first term
%
if isa(bc.term1, 'function_handle')
  [rV1, rV2] = self.callValueFuncOnEdge(bc, femEdges, bc.term1);
  nrRV1 = size(rV1,1); nrRV2 = size(rV2,1);
  % The number of constrained dofs at both ends of the edge can't
  % be larger than the total number of dofs
  if(nrRV1 > self.N)
    eId = 'pde:pde2DBCImpl:t1TooLarge';
    err = message(eId, nrRV1, self.N, edgeID);
    throwAsCaller(MException(eId, err.getString()));
  end
  if(nrRV2 > self.N)
    eId = 'pde:pde2DBCImpl:t1TooLarge';
    err = message(eId, nrRV2, self.N, edgeID);
    throwAsCaller(MException(eId, err.getString()));
  end
  % Make sure the lengths of term1 and term2 are compatible at both
  % ends of the edge.
  checkTermLengths(nrRV1, lenRr1, bc.appRegionID);
  checkTermLengths(nrRV2, lenRr1, bc.appRegionID);
  r(rr1,femEdges) = rV1;
  r(rr1,femEdges+ne) = rV2;
elseif isnumeric(bc.term1)
  ri = zeros(length(rr1), 1);
  if ~isempty(bc.term1)
    lenT1 = length(bc.term1);
    if(lenT1 > self.N)
      eId = 'pde:pde2DBCImpl:t1TooLarge';
      err = message(eId, lenT1, self.N, bc.appRegionID);
      throwAsCaller(MException(eId, err.getString()));
    end
    checkTermLengths(lenT1, lenRr1, bc.appRegionID);
    ri = zeros(length(rr1), 1);
    ri(:) = bc.term1(:);
  end
  rAllEdges = repmat(ri, 1, numFemEdges);
  r(rr1,femEdges) = rAllEdges;
  r(rr1,femEdges+ne) = rAllEdges;
else
  error(message('pde:pde2DBCImpl:term1NotNum', edgeID));
end

end

function hN = getHInd(N, dofList, edgeID)
% get indices into  an N*N h-vector for the list of dofs
lenDofs = length(dofList);
hN = zeros(lenDofs, 1);
for k=1:lenDofs
  kk = dofList(k);
  if(kk < 1 || kk > N || ~isreal(kk) || mod(kk,1) ~= 0)
    kkStr = num2str(kk, 16);
    error(message('pde:pde2DBCImpl:eqnNumInvalidValue', kkStr, edgeID, N));
  end
  hN(k) = (kk-1)*N + kk;
end
end

function checkTermLengths(lenT1, lenT2, edgeID)
% The lengths of term1 and term2 must be the same
if(lenT1~=lenT2)
  eId = 'pde:pde2DBCImpl:valUnequalLengths';
  err = message(eId, lenT1, lenT1, lenT2, edgeID);
  throwAsCaller(MException(eId, err.getString()));
end
end

