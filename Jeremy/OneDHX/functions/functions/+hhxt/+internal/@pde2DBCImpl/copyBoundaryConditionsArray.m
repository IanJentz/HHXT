function self = copyBoundaryConditionsArray(self)
% copy pdeBoundaryConditions to local data structure
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014-2015 The MathWorks, Inc.
import pde.internal.pdeBCImpl;

if(~isempty(self.myPde))
  % Use a map to collect all BCs for an edge with a given integer ID
  self.edgeBoundaryConditions = containers.Map('KeyType','uint32','ValueType','any');
  for i=1:length(self.myPde.BoundaryConditions)
    bci = self.myPde.BoundaryConditions(i);
    % A case with no parameters is a zero-Neumann BC 
    % This is the default for FEM so we just skip it
    if isempty(bci.Type)
      continue;
    end
    if(strcmpi(bci.Type,'neumann'))
      val.bcType = self.neumannBC;
      val.term1 = bci.g;
      val.term2 = bci.q;
    elseif(strcmpi(bci.Type, 'dirichlet'))
      val.bcType = self.dirichletBC;
      val.term1 = bci.r;
      val.term2 = bci.h;
    elseif(strcmpi(bci.Type, 'value'))
      val.bcType = self.valueBC;
      val.term1 = bci.u;
      val.term2 = bci.EquationIndex;
    else
      error(message('pde:pde2DBCImpl:invalidBCType', bci.Type));
    end
    val.isVectorized = strcmpi(bci.Vectorized, 'on');
        
    % Add all boundary conditions to map with key as edge id
	% This map will be accessed in the routines that calculate the BC coefficient matrices    
    for j=1:numel(bci.RegionID)
      val.appRegionID = bci.RegionID(j);
      jj = bci.RegionID(j);
      if strcmp(bci.RegionType, 'Edge')   
        if(self.edgeBoundaryConditions.isKey(jj))
          self.edgeBoundaryConditions(jj) = [self.edgeBoundaryConditions(jj) val];
        else
          self.edgeBoundaryConditions(jj) = val;
        end
      end
    end
    
    
  end
end
end

