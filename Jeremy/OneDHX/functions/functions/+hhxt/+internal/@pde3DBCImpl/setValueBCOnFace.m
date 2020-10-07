function setValueBCOnFace( self, bci, dirBoundaryConditions )
%setValueBCOnFace Calculate r and h entries for a value type BC
% This undocumented function may be changed or removed in a future release.
% This function handles the case when a 'u' and/or 'EquationIndex' parameter
% has been included in a pdeBoundaryConditions entity.

%       Copyright 2014-2016 The MathWorks, Inc.

import pde.internal.*;

faceID = bci.appRegionID;
numNodes = length(bci.nodes);

% Both the first and second terms are optional. But to get to this point,
% at least one must have been included to determine the type.
nodesCell = num2cell(bci.nodes);
doesNodeHaveBC = dirBoundaryConditions.isKey(nodesCell);

% Following flags are used to check if any node on this FE face
% or any DoF at a node still needs BC assignment. If not, do
% nothing and return to the caller.

nodesNeedBCAssignment = false;
dofsNeedBCAssignment = false;
if (~all(doesNodeHaveBC))
    nodesNeedBCAssignment = true;
end

if (nodesNeedBCAssignment == false)
    % check for BC at all DoFs only when all nodes already have some BC
    for n = cell2mat(nodesCell(doesNodeHaveBC))'
        hr = dirBoundaryConditions(n);
        previousBCAssignment = hr.r ~=0;
        for j = 1:numel(hr.r)
            if (~previousBCAssignment(j)) % If there was no previous assignment
                dofsNeedBCAssignment = true;
                continue
            end
        end
    end
end

if ~(nodesNeedBCAssignment || dofsNeedBCAssignment)
    return
end

isFuncU = isa(bci.term1, 'function_handle');

% Check to see if a second term has been included.
if(~isempty(bci.term2))
    % we have a second term
    if ~isnumeric(bci.term2)
        error(message('pde:pde2DBCImpl:term2NotNum', faceID));
    end
    t2 = bci.term2;
    rr1 = t2(:);
    lenRr1 = length(rr1);
    ri = zeros(self.N,1);
    rnz = ri~=0;
    if(lenRr1 > self.N)
        eId = 'pde:pde2DBCImpl:t2TooLarge';
        err = message(eId, lenRr1, self.N, faceID);
        throwAsCaller(MException(eId, err.getString()));
    end
    hMat = getH(self.N, t2, faceID);
else
    % no second term, default h-matrix (all dofs)
    rr1 = 1:self.N;
    lenRr1 = self.N;
    hMat = eye(self.N);
end
%
% Now handle the first term
%

if isFuncU
    %    if ~all(doesNodeHaveBC)
    faceNormals = self.calcFaceNormalsAtNodes(bci.elemFaces, bci.nodes);
    uVec = self.callValueFuncOnFace(bci, bci.nodes, bci.term1, faceNormals);
    nr = size(uVec, 1);
    if(nr > self.N)
        eId = 'pde:pde2DBCImpl:t1TooLarge';
        err = message(eId, nr, self.N, faceID);
        throwAsCaller(MException(eId, err.getString()));
    end
    checkTermLengths(nr, lenRr1, faceID);
    rMat = zeros(self.N, numNodes);
    rMat(rr1,:) = uVec;
    %    end
    
elseif isnumeric(bci.term1)
    rMat = zeros(self.N, 1);
    if ~isempty(bci.term1)
        lenT1 = length(bci.term1);
        if(lenT1 > self.N)
            eId = 'pde:pde2DBCImpl:t1TooLarge';
            err = message(eId, lenT1, self.N, faceID);
            throwAsCaller(MException(eId, err.getString()));
        end
        checkTermLengths(lenT1, lenRr1, faceID);
        rMat(rr1) = bci.term1(:);
        ri = rMat;
        rnz = ri~=0;
    end
else
    error(message('pde:pde2DBCImpl:term1NotNum', faceID));
end

hrZ = HRNode(self.N);
hi = hMat;
hnz = hi~=0;


for i=1:numNodes
    nodeID = bci.nodes(i);
    if(doesNodeHaveBC(i))
        hr = dirBoundaryConditions(nodeID);
    else
        hr = hrZ;
    end
    if(isFuncU)
        ri = rMat(:,i);
        rnz = ri~=0;
    end
    hr.nodeID = nodeID;
    hr.r(rnz) = ri(rnz);
    hr.h(hnz) = hi(hnz);
    dirBoundaryConditions(nodeID) = hr;
end
end


function hN = getH(N, dofList, edgeID)
% get indices into  an N*N h-vector for the list of dofs
lenDofs = length(dofList);
hN = zeros(N);
for k=1:lenDofs
    kk = dofList(k);
    if(kk < 1 || kk > N || ~isreal(kk) || mod(kk,1) ~= 0)
        kkStr = num2str(kk, 16);
        error(message('pde:pde2DBCImpl:eqnNumInvalidValue', kkStr, edgeID, N));
    end
    hN(kk,kk) = 1;
end
end

function checkTermLengths(lenT1, lenT2, edgeID)
if(lenT1~=lenT2)
    eId = 'pde:pde2DBCImpl:valUnequalLengths';
    err = message(eId, lenT1, lenT1, lenT2, edgeID);
    throwAsCaller(MException(eId, err.getString()));
end
end