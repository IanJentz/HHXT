classdef HHXTModel < pde.PDEModel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
      properties
        % EquationCoefficients - Equation coefficients within the domain
        %    An object that contains the equation coefficient assignments within the
        %    geometric domain. Equation coefficients are specified and assigned to
        %    the geometry using the specifyCoefficients function. If multiple assignments
        %    are made, they are all recorded within this object.
        %
        %    See also pde.CoefficientAssignmentRecords/findCoefficients,
        %             pde.CoefficientAssignmentRecords/globalCoefficients
        AddtEquationCoefficients;
        
        MaterialProperties;
        
        MatrixOptions;
        
        HHXTSolverOptions;
        
        MaterialMap
        
        AnalysisType;
        
        WorkingDirectory;
        
        Thickness;
        
        AdvectTermDemotion;
        
        PermeabilityUpdateCount;
        
        Runtime;
                
        % StreamBoundaryConditions - Boundary conditions applied to the
        % streams of the HHXT geometry
        %    An object that contains the boundary condition assignments on the
        %    boundary of the geometric. Boundary conditions are specified and assigned to
        %    the geometry using the applyBoundaryCondition function. If multiple assignments
        %    are made, they are all recorded within this object.
        %
        %    See also pde.BoundaryConditionRecords/findBoundaryConditions
        StreamBoundaryConditions;
        
      end
    
    methods(Static)
        
%         sol = solveHHXTpde(self, varargin) % edited from pde toolbox 
    
    end
    
    methods
        function obj = HHXTModel(numEqns)
            if(nargin < 1)
                numEqns = 1;
            end
            validateattributes(numEqns,{'numeric'},...
                {'scalar','integer','positive'});
            obj@pde.PDEModel(numEqns)
            
            obj.PermeabilityUpdateCount = 0;
            obj.Thickness = 1;
            
        end
        function obj = setThickness(obj,H)
            obj.Thickness = H;
        end
   % Method declaration
        coef = specifyHHXTCoefficients(self, varargin)
        sol = solveHHXTpde(self, varargin) % edited from pde toolbox 
        u = solveStationaryHHXT(self,coefstruct,coefstructb,varargin) % edited from pde.EquationModel
        u = solveStationaryNonlinearHHXT(self, coefstruct, coefstructb, u0) % edited from pde.EquationModel
        [u,dudt,dudt2,tlist] = solveTimeDependentHHXT(self, coefstruct, coefstructb, u0, ut0, tlist, tsecondOrder) % edited from pde.EquationModel
        self = updatePerm(self)
        self = modelQuadPoints(self)
        self = updateQuadPointVals(self,u,p,e,t)
        
    end
    
    
end

