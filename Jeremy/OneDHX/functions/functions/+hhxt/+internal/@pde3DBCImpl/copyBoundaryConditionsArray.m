function self = copyBoundaryConditionsArray(self,bc,isR2016b)
% copy pdeBoundaryConditions to local data structure
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014-2015 The MathWorks, Inc.

import pde.internal.pdeBCImpl;

if(~isempty(self.myPde))
    self.faceBoundaryConditions = struct('bcType',{},'term1',{},'term2',{}, ...
        'isVectorized',{},'appRegionID',{}, 'faceID',{});
    if isempty(bc)
        self.faceBoundaryConditions = [];
    elseif isR2016b
        for i=1:length(bc)
            bci = bc(i);
            val = struct();
            if(strcmpi(bci.Type,'neumann'))
                val.bcType = pdeBCImpl.neumannBC;
                val.term1 = bci.g;
                val.term2 = bci.q;
            elseif(strcmpi(bci.Type,'dirichlet'))
                val.bcType = pdeBCImpl.dirichletBC;
                val.term1 = bci.r;
                val.term2 = bci.h;
            elseif(strcmpi(bci.Type,'value'))
                val.bcType = pdeBCImpl.valueBC;
                val.term1 = bci.u;
                val.term2 = bci.EquationIndex;
            elseif(strcmpi(bci.Type,'mixedvalue'))
                val(1).bcType = pdeBCImpl.valueBC;
                val(1).term1 = bci.u;
                val(1).term2 = bci.EquationIndex;
                val(2).bcType = pdeBCImpl.neumannBC;
                val(2).term1 = bci.g;
                val(2).term2 = bci.q;
            elseif(strcmpi(bci.Type,'mixeddirichlet'))
                val(1).bcType = pdeBCImpl.dirichletBC;
                val(1).term1 = bci.r;
                val(1).term2 = bci.h;
                val(2).bcType = pdeBCImpl.neumannBC;
                val(2).term1 = bci.g;
                val(2).term2 = bci.q;
            else
                error(message('pde:pde2DBCImpl:invalidBCType', bci.Type));
            end
            for j = 1:numel(val)
                val(j).isVectorized = strcmpi(bci.Vectorized, 'on');
                val(j).appRegionID = i;
                val(j).faceID = i;
                self.faceBoundaryConditions(end+1) = val(j);
            end
        end
    else
        for i=1:length(bc)
            bci = bc(i);
            % A case with no parameters is a zero-Neumann BC
            % This is the default for FEM so we just skip it
            if isempty(bci.Type)
                continue;
            end
            if(strcmpi(bci.Type,'neumann'))
                val.bcType = pdeBCImpl.neumannBC;
                val.term1 = bci.g;
                val.term2 = bci.q;
            elseif(strcmpi(bci.Type, 'dirichlet')||strcmpi(bci.Type, 'mixeddirichlet'))
                val.bcType = pdeBCImpl.dirichletBC;
                val.term1 = bci.r;
                val.term2 = bci.h;
            elseif(strcmpi(bci.Type, 'value')||strcmpi(bci.Type, 'mixedvalue'))   
                val.bcType = pdeBCImpl.valueBC;
                val.term1 = bci.u;
                val.term2 = bci.EquationIndex;
            else
                error(message('pde:pde2DBCImpl:invalidBCType', bci.Type));
            end
            val.isVectorized = strcmpi(bci.Vectorized, 'on');
            
            
            % A pdeBoundaryCondition can have multiple faces in its application region
            % Separate them out here.
            
            for j=1:numel(bci.RegionID)
                val.appRegionID = bci.RegionID(j);
                if strcmp(bci.RegionType, 'Face')
                    jj = bci.RegionID(j);
                    val.faceID = jj;
                    self.faceBoundaryConditions(end+1) = val;
                end
            end
            
        end
        
    end
end


end

