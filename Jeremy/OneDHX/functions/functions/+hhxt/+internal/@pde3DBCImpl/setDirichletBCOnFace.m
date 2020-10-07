function setDirichletBCOnFace( self, bci, dirBoundaryConditions )
%setDirichletBCOnFace Calculate r and h entries for a Dirichlet BC
% This undocumented function may be changed or removed in a future release.
% This function handles the case when a 'r' and/or 'h' parameter has been
% included in a pdeBoundaryConditions entity.

%       Copyright 2014-2016 The MathWorks, Inc.

import pde.internal.*;

faceID = bci.appRegionID;
numNodes = length(bci.nodes);

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

% Handle the first term (r vector)
isFuncR = isa(bci.r, 'function_handle');
isFuncH = isa(bci.h, 'function_handle');
if(isFuncR || isFuncH)
    faceNormals = self.calcFaceNormalsAtNodes(bci.elemFaces, bci.nodes);
end
if(~isempty(bci.term1))
    if isFuncR
        rMat = self.callDirichletFuncOnFace(bci, bci.nodes, bci.r, faceNormals, 1);
    else
        ri = bci.term1(:);
        if(length(ri) ~= self.N)
            error(message('pde:pde2DBCImpl:invalidLenR', faceID, self.N));
        end
        rnz = ri~=0;
    end
else
    ri = zeros(self.N,1);
    rnz = ri~=0;
end

% Handle the second term (h matrix)
if(isempty(bci.term2))
    hi = eye(self.N);
    hnz = hi~=0;
else
    if isFuncH
        hMat = self.callDirichletFuncOnFace(bci, bci.nodes, bci.h, faceNormals, 2);
        hMat = reshape(hMat, self.N, self.N, []);
    else
        hi = bci.term2(:);
        if(length(hi) ~= self.N^2)
            error(message('pde:pde2DBCImpl:invalidLenH', faceID, self.N^2));
        end
        hnz = hi~=0;
    end
end

hrZ= HRNode(self.N);

for i=1:numNodes
    nodeID = bci.nodes(i);
    if(doesNodeHaveBC(i))
        hr = dirBoundaryConditions(nodeID);
    else
        hr = hrZ;
    end
    if(isFuncH)
        hi = hMat(:,:,i);
        hnz = hi~=0;
    end
    if(isFuncR)
        ri = rMat(:,i);
        rnz = ri~=0;
    end
    hr.r(rnz) = ri(rnz);
    hr.h(hnz) = hi(hnz);
    dirBoundaryConditions(nodeID) = hr;
end
end


