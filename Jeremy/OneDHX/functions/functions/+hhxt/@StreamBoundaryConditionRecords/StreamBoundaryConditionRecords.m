classdef (Sealed) StreamBoundaryConditionRecords  < handle & matlab.mixin.internal.CompactDisplay
% BoundaryConditionRecords  Assignments of the boundary conditions
%    The BoundaryConditionRecords holds a record of all the assignments of 
%    boundary conditions across the geometric boundary. This may be a single
%    assignment that represents the same boundary conditions throughout the 
%    entire boundary or multiple assignments where the boundary conditions 
%    are different across various boundaries. 
%
% BoundaryConditionRecords methods:
%    findBoundaryConditions  - Find boundary condition assignment for a
%    boundary
%
% BoundaryConditionRecords properties:
%    BoundaryConditionAssignments - Vector of boundary condition assignments
%
% See also hhxt.HHXTModel, hhxt.HHXTModel/applyStreamBoundaryCondition

% Copyright 2020 Ian Jentz

properties
    
% BoundaryConditionAssignments - Vector of boundary condition assignments
%    A vector containing the boundary condition assignments to the boundaries.
%    Each entry in the vector is a pde.BoundaryCondition object. 
%    This object defines the boundary condition of the PDE. Coefficients and 
%    host boundaries are defined using the applyBoundaryCondition
%    method of the PDEModel.
    BoundaryConditionAssignments;
end

methods
   % Method declaration
   bc = findBoundaryConditions(self, varargin)   
end

methods (Hidden = true, Access = private)
    tf = boundaryConditionsSpecified(self, varargin)  
end

methods(Hidden=true)
    function obj=StreamBoundaryConditionRecords(pdem)
       obj.BoundaryConditionAssignments = [];
       obj.ParentPdemodel = pdem;
    end 
end

methods(Hidden=true, Access={?hhxt.HHXTModel, ?hhxt.HHXTResults, ?hhxt.StreamBoundaryCondition})
    function bcArray = packBoundaryConditions(self)
        if isempty(self)
            bcArray = [];
            return
        end
        if self.ParentPdemodel.IsTwoD
            nbc = self.ParentPdemodel.Geometry.NumEdges;
            region = 'edge';
        else
            nbc = self.ParentPdemodel.Geometry.NumFaces;
            region = 'face';
        end
        tf = boundaryConditionsSpecified(self, region, 1:nbc);
        allbc = 1:nbc;
        nbctrue = allbc(tf);
        bc = hhxt.StreamBoundaryCondition(self,'HeatFlux',0,region,1);
        bcArray = repmat(bc,nbc,1);
        for i = 1:length(nbctrue)
            bc = self.findBoundaryConditions(region,nbctrue(i));
            bcArray(nbctrue(i)) = bc;
        end
        
        % compile the stream BCs into matlab pde BCs
        N = self.ParentPdemodel.PDESystemSize;
        
        nbcs = size(self.BoundaryConditionAssignments,2);
        rbcs = ones(1,nbcs);
        for i = 1:nbcs
           rbcs(i) = self.BoundaryConditionAssignments(i).RegionID; 
        end
        
        bcsfreg = cell(1,length(allbc));
        for i = allbc
            bcsfreg{i} = find(rbcs==i);
        end
        
        for j = unique(rbcs)
            rm = zeros(N,1);
            hm = zeros(N,1);
            gm = zeros(N,1);
            nDirichlet = 0;
            nNeumann = 0;
            for i = bcsfreg{j}
                ntmp = size(self.BoundaryConditionAssignments(i).ThermalTypes,2);
                nhyd = size(self.BoundaryConditionAssignments(i).HydraulicTypes,2);
                str = self.BoundaryConditionAssignments(i).StreamID;
                area = self.BoundaryConditionAssignments(i).BCarea;
                dirmult = 1;
                if strcmp(self.BoundaryConditionAssignments(i).FlowDir,'outlet')
                    dirmult = -1;
                end
                for nt = 1:ntmp
                    val = self.BoundaryConditionAssignments(i).ThermalVals{nt};
                    switch self.BoundaryConditionAssignments(i).ThermalTypes{nt}
                        case 'Temperature'
                            nDirichlet = nDirichlet + 1;
                            hm(2*str+1) = 1;
                            rm(2*str+1) = rm(2*str+1) + val; 
                        case 'HeatFlux'
                            nNeumann = nNeumann + 1;
                            gm(2*str+1) =  + val;
                        case 'HeatFlow'
                            nNeumann = nNeumann + 1;
                            gm(2*str+1) =  + dirmult*val/area;
                    end
                end
                for nh = 1:nhyd
                    val = self.BoundaryConditionAssignments(i).HydraulicVals{nh};
                    switch self.BoundaryConditionAssignments(i).HydraulicTypes{nh}
                        case 'Pressure'
                            nDirichlet = nDirichlet + 1;
                            hm(2*str) = 1;
                            rm(2*str) = rm(2*str) + val; 
                        case 'MassFlux'
                            nNeumann = nNeumann + 1;
                            gm(2*str) =  + dirmult*val;
                        case 'MassFlow'
                            nNeumann = nNeumann + 1;
                            gm(2*str) =  + dirmult*val/area;
                    end
                end
            end
            hm = diag(hm);
            if ( (nDirichlet ~= 0) && (nNeumann ~= 0))
                formulation = 'mixed';
                applyBoundaryCondition(self.ParentPdemodel,'mixed',region,j,'g',gm,'h',hm,'r',rm);
            elseif nDirichlet ~=0
                formulation = 'dirichlet';
                applyBoundaryCondition(self.ParentPdemodel,'dirichlet',region,j,'h',hm,'r',rm);
            else
                formulation = 'neumann';
                applyBoundaryCondition(self.ParentPdemodel,'neumann',region,j,'g',gm);
            end
            for k = bcsfreg{j}
                self.BoundaryConditionAssignments(k).BCType = formulation;
            end
            
            
        end
        
        
    end
end
  
methods(Hidden=true, Access={?pde.BoundaryCondition,?pde.EquationModel})
    function delistBoundaryConditions(self, bctodelete)
        numbc = numel(self.StreamBoundaryConditionAssignments);
        for i = 1:numbc
            thisbc = self.StreamBoundaryConditionAssignments(i);
            if thisbc == bctodelete
                self.StreamBoundaryConditionAssignments(i) = [];
                break
            end
        end  
        numbc = numel(self.StreamBoundaryConditionAssignments);
        if numbc == 0 && ~isempty(self.ParentPdemodel) && isvalid(self.ParentPdemodel)
           self.ParentPdemodel.delistBoundaryConditions(); 
        end
    end 
end

methods(Hidden=true)
    function delete(self)
        numic = numel(self.BoundaryConditionAssignments);
        for i = 1:numic
            if isvalid(self.BoundaryConditionAssignments(i))           
                delete(self.BoundaryConditionAssignments(i));
            end
        end   
               
        if isvalid(self.ParentPdemodel)
            self.ParentPdemodel.delistBoundaryConditions();
        end        
    end    
end

properties (Hidden = true, SetAccess='private')
    ParentPdemodel;
end  
  
end
