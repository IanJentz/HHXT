classdef(Hidden) pde2DBCImpl < pde.internal.pdeBCImpl
  %pde2DBCImpl Implement a pdebound function for 2D boundary conditions
  %
  % This undocumented function may be changed or removed in a future release.
  
  %       Copyright 2014 The MathWorks, Inc.
  
  properties(Access=private)
    edgeBoundaryConditions;
    points, numPoints, edges, numEdges, time, uN, uxN, uyN, gotTime, gotU;
    hDefault;
  end
  
  methods
    function obj = pde2DBCImpl(myPde)
      obj@pde.internal.pdeBCImpl(myPde);
      obj.copyBoundaryConditionsArray();
    end
    function fh = getBoundaryFunc(self)
      % Implements a function handle to a "pdebound" function for
      % boundary conditions on 2D geometry. Used by assempde, hyperbolic,
      % etc.
      fh = @self.getBCMatrices;
    end
  end % methods
  
  methods(Access=private)
    % Class methods defined outside of this file.
    [q,g,h,r] = getBCMatrices(self,p,e,u,time);
    copyBoundaryConditionsArray(self);
    [h, r] = setValueBCOnEdge(self, bc, femEdges, h, r);
    [q, g] = setNeumannBCOnEdge(self, bc, femEdges, q, g);
    [h, r] = setDirichletBCOnEdge(self, bc, femEdges, h, r);
    [rhV1, rhV2] = callDirichletFuncOnEdge(self, bc, femEdges, f, numRows);
    qg = callNeumannFuncOnEdge(self, bc, femEdges, f, numRows);
    [v1, v2] = callValueFuncOnEdge(self, bc, femEdges, f);
    [nx, ny] = calcEdgeNormals(self, femEdges);
  end % methods
  
end

