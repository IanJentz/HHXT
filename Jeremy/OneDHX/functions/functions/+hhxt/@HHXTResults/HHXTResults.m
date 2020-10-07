classdef HHXTResults < pde.StationaryResults
    % pde.StationaryResults PDE solution and its derived quantities
    %   A StationaryResults object provides a convenient representation of
    %   results data for stationary PDE analysis, together with
    %   interpolation and gradient evaluation functions. Create a
    %   StationaryResults object using the createPDEResults function.
    %
    % StationaryResults methods:
    %   interpolateSolution  - Interpolate solution at specified spatial
    %                          locations
    %   evaluateGradient     - Evaluate gradients of solution at specified
    %                          spatial locations
    %   evaluateCGradient    - Evaluate tensor product of c-coefficient and
    %                          gradients of the PDE solution
    %
    % StationaryResults properties:
    %   NodalSolution        - Solution to PDE at nodal locations
    %   XGradients           - Spatial gradient of solution along x-direction.
    %   YGradients           - Spatial gradient of solution along y-direction.
    %   ZGradients           - Spatial gradient of solution along z-direction.
    %   Mesh                 - Discretization of the domain
    %
    %    See also createPDEResults
    
    % Copyright 2015-2016 The MathWorks, Inc.
    
    
    properties(SetAccess = protected)
        
        DarcyFlux;
        XDarcyFlux;
        YDarcyFlux;
        ZDarcyFlux;
        ChannelVelocity;
        XChannelVelocity;
        YChannelVelocity;
        ZChannelVelocity;
        Reynolds;
        XReynolds;
        YReynolds;
        ZReynolds;
        StreamVolHeating;
        StreamVolResistance;
        
        ElemHeating;
        
        Tables

    end
    
    methods
        function obj = HHXTResults(varargin)
            
            pdem = varargin{1};
            u = varargin{2};
            u(isnan(u)) = 0; % maybe add a NaN check here? 
            varargin{2} = u;
            
            calcVars = true;
            if nargin == 3
                calcVars = varargin{3};
                varargin(3) = [];
            end
            
            obj@pde.StationaryResults(varargin{:});
            if nargin == 0
                return
            end
            narginchk(2,3);
            if (obj.IsTimeEig)
                error(message('pde:PDEResults:notStationaryResults'));
            end
                       
            ureshaped = pde.PDEResults.reshapePDESolution(u, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            obj.NodalSolution = ureshaped;
            if ~isempty(pdem.EquationCoefficients)
                obj.Coefficients = pdem.EquationCoefficients.packCoefficients;
            else
                obj.Coefficients = [];
            end
            
            obj.XGradients = [];
            obj.YGradients = [];
            obj.ZGradients = [];
            obj.InterpolantdUdx = [];
            obj.InterpolantdUdy = [];
            obj.InterpolantdUdz = [];
            if isempty(obj.NodalSolution)
                return;
            end            
            
            obj.Interpolant   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                ureshaped, pdem.PDESystemSize, obj.NumTimeEig);
            
            % calculate gradients, the HHXT way
            % this allows for internal boundaries to exist
            % creates obj.XGradients, obj.YGradients, obj.ZGradients
            obj = calcGradients(obj,pdem);
              
            obj.InterpolantdUdx   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.XGradients, pdem.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdy   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.YGradients, pdem.PDESystemSize, obj.NumTimeEig);
            
            if ~obj.IsTwoD
                obj.InterpolantdUdz   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZGradients, pdem.PDESystemSize, obj.NumTimeEig);
            end
            
            if calcVars == true
            
            % calculate stream variables, Darcy Velocity, channel velocity,
            % reynolds number, thermal restance, and heating
            obj = calcVariables(obj,pdem);
            
            nstream = (pdem.PDESystemSize-1)/2;
            obj.InterpolantDarcyFlux   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.DarcyFlux, nstream, obj.NumTimeEig);
            obj.InterpolantXDarcyFlux   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.XDarcyFlux, nstream, obj.NumTimeEig);
            obj.InterpolantYDarcyFlux   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.YDarcyFlux, nstream, obj.NumTimeEig);
            
            obj.InterpolantChannelVelocity   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.ChannelVelocity, nstream, obj.NumTimeEig);
            obj.InterpolantXChannelVelocity   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.XChannelVelocity, nstream, obj.NumTimeEig);
            obj.InterpolantYChannelVelocity   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.YChannelVelocity, nstream, obj.NumTimeEig);
            
            obj.InterpolantReynolds   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.Reynolds, nstream, obj.NumTimeEig);
            obj.InterpolantXReynolds   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.XReynolds, nstream, obj.NumTimeEig);
            obj.InterpolantYReynolds   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.YReynolds, nstream, obj.NumTimeEig);
            
            if ~obj.IsTwoD
                obj.InterpolantZDarcyFlux   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZDarcyFlux, nstream, obj.NumTimeEig);
                obj.InterpolantZChannelVelocity   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZChannelVelocity, nstream, obj.NumTimeEig);
                obj.InterpolantZReynolds   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZReynolds, nstream, obj.NumTimeEig);
            end
            
            obj.InterpolantStreamVolHeating   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.StreamVolHeating, nstream, obj.NumTimeEig);
            
            obj.InterpolantStreamVolResistance   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.StreamVolResistance, nstream, obj.NumTimeEig);
            
            end
            
            obj = calcBoundaries(obj,pdem);
            
        end
        
        % Methods declaration
        obj = calcGradients(obj,hhxtm)
        obj = calcVariables(obj,hhxtm)
        obj = calcBoundaries(obj,hhxtm)
        
    end
    
    properties (Hidden = true, Transient=true)
        InterpolantDarcyFlux;
        InterpolantXDarcyFlux;
        InterpolantYDarcyFlux;
        InterpolantZDarcyFlux
        InterpolantChannelVelocity;
        InterpolantXChannelVelocity;
        InterpolantYChannelVelocity;
        InterpolantZChannelVelocity;
        InterpolantReynolds;
        InterpolantXReynolds;
        InterpolantYReynolds;
        InterpolantZReynolds;
        InterpolantStreamVolHeating;
        InterpolantStreamVolResistance;
    end
    
      
    
end


