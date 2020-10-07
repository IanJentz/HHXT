classdef(Hidden) pde3DBCImpl < pde.internal.pdeBCImpl
  %pde3DBCImpl Implement a pdebound function for 3D boundary conditions
  %
  % This undocumented function may be changed or removed in a future release.
  
  %       Copyright 2014-2015 The MathWorks, Inc.
  
  properties(Access=private)
    faceBoundaryConditions;
    points, numPoints, edges, numEdges, time, uN, uxN, uyN, gotTime, gotU;
    numDirichletNodes;
  end
  
  methods
    function obj = pde3DBCImpl(myPde,bc,isR2016b)
      obj@pde.internal.pdeBCImpl(myPde);
      obj.copyBoundaryConditionsArray(bc,isR2016b);
      obj.numDirichletNodes = 0;
    end
    [Q,G,R,H] = getBCMatrices(self,p,geomMeshAssoc,u,time);
  end % methods
  
  methods(Access=private)
    % Class methods defined outside of this file.
    copyBoundaryConditionsArray(self,bc,isR2016b);
    setValueBCOnFace(self, bc, dirBoundaryConditions);
    [Qi, Gi] = setNeumannBCOnFace(self, bc);
    setDirichletBCOnFace(self, bci, dirBoundaryConditions);
    rh = callDirichletFuncOnFace(self, bc, nodes, func, faceNormals, t12);
    qg = callNeumannFuncOnFace(self, bc, xyz, sPts, func, faceNormVec, uPoints, t12);
    rMat = callValueFuncOnFace(self, bc, nodes, func, faceNormals);
    faceNormals = calcFaceNormalsAtNodes(self, faces, nodes);
  end % methods
  
  methods(Static, Access=private)
    function nodes=getNodesForFaces(faces)
      nodes = unique(faces(:));
    end
    appRegion = applicationRegion(locations, normals);
  end % methods
  
end

