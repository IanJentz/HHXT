function faceNormVecAtNodes = calcFaceNormalsAtNodes(self, faces, nodes)
%calcFaceNormalsAtNodes Calculate normal to a face at each of its nodes
% This undocumented function may be changed or removed in a future release.
% Given a set of faces and face nodes, this function calculates the face 
% normal at each node. 

%       Copyright 2014 The MathWorks, Inc.

import pde.internal.elements.*;


[numFaceNodes, numFaces] = size(faces);

% Get the shape function derivatives wrt the natural coordinates
% for the triangle shape that represents a face of a tetrahedron.
if(numFaceNodes == 3)
  dsfCalc = @dShapeTri3;
  rs = [0 1 0; 0 0 1]; % natural coordinates at the 3 nodes
  numNodesToEval = 1; % for flat elements, compute normal once
elseif(numFaceNodes == 6)
  dsfCalc = @dShapeTri6;
  rs = [0 1 0 .5 .5 0; 0 0 1 0 .5 .5]; % natural coordinates at the 6 nodes
  numNodesToEval = 6;
else
  error(message('pde:pde3DBCImpl:illegalFaceDefn', numFaceNodes));
end
% Collect the shape function derivatives for all face nodes 
dsNodes = zeros(numFaceNodes, 2, numFaceNodes);
for n=1:numFaceNodes
  dsNodes(:,:,n) = dsfCalc(rs(:,n));
end

numNodes = length(nodes);
numFacesAtNode = zeros(numNodes, 1);
faceNormVecAtNodes = zeros(3, numNodes);

nodeIndexVec = zeros(self.numPoints, 1);
nodeIndexVec(nodes) = 1:numNodes;

xyzAllFaceNodes = self.points(:, faces(:));
xyzAllFaceNodes = reshape(xyzAllFaceNodes, 3, numFaceNodes, numFaces);
faceNormVec = zeros(3, numFaceNodes);
for i=1:numFaces
  xyz = xyzAllFaceNodes(:,:,i);
  for j=1:numNodesToEval
    % Compute the two vectors dX/dr and dX/ds
    dxyzDrs = xyz*dsNodes(:,:,j);
    % Their cross product is the normal at the point
    v3 = cross3(dxyzDrs(:,1),dxyzDrs(:,2));
    faceNormVec(:,j) = v3/norm(v3);
  end
  if(numNodesToEval==1)
    faceNormVec = repmat(faceNormVec(:,1), 1, numFaceNodes);
  end
  for j=1:numFaceNodes
    nj = nodeIndexVec(faces(j,i));
    numFacesAtNode(nj) = numFacesAtNode(nj) + 1;
    faceNormVecAtNodes(:,nj) = faceNormVecAtNodes(:,nj) + faceNormVec(:,j);
  end
end

% The normal at a node is the average of the normals of all the facets
% connected to that node.
for i=1:numNodes
  faceNormVecAtNodes(:,i) = faceNormVecAtNodes(:,i)/numFacesAtNode(i);
end

end

% This highly specialized version of the cross function with no error
% checking is used to improve performance.
function c=cross3(a,b) 
c=[a(2)*b(3)-a(3)*b(2); 
   a(3)*b(1)-a(1)*b(3); 
   a(1)*b(2)-a(2)*b(1)];
end

