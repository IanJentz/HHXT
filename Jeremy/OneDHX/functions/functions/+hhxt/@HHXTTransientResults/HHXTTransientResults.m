classdef HHXTTransientResults < pde.PDEResults
    % pde.TimeDependentResults PDE solution and its derived quantities
    %   A TimeDependentResults object provides a convenient representation
    %   of results data for a time-dependent PDE analysis, together with
    %   interpolation and gradient evaluation functions. Create a
    %   TimeDependentResults object using the createPDEResults function.
    %
    % TimeDependentResults methods:
    %   interpolateSolution  - Interpolate solution at specified spatial
    %                          locations
    %   evaluateGradient     - Evaluate gradients of solution at specified
    %                          spatial locations
    %   evaluateCGradient    - Evaluate tensor product of c-coefficient and
    %                          gradients of the PDE solution    
    %
    % TimeDependentResults properties:
    %   NodalSolution        - Solution to PDE at nodal locations
    %   SolutionTimes        - Times at which the NodalSolution was computed
    %   XGradients           - Spatial gradient of solution along x-direction
    %   YGradients           - Spatial gradient of solution along y-direction
    %   ZGradients           - Spatial gradient of solution along z-direction
    %   Mesh                 - Discretization of the domain
    %
    %    See also createPDEResults
    
    % Copyright 2015-2016 The MathWorks, Inc.
    
    
    properties(SetAccess = protected)
        % NodalSolution - Solution to PDE at nodal locations
        % Solution array which can be conveniently indexed into to extract
        % a sub-array of interest. The shape of NodalSolution depends on
        % the type of PDE and solver settings. It will be a:
        %
        %   column vector - for a single PDE with no time dependency
        %   matrix        - for a single hyperbolic or parabolic problem,
        %                   or a system of elliptic problems, or a single
        %                   eigenvalue problem
        %   3-D array     - for a system of hyperbolic, parabolic, or
        %                   eigenvalue problems
        %
        % The first array dimension of NodalSolution represents node index.
        % The second array dimension represents the time-step or
        % eigenvector index for a single PDE, or the equation index for a
        % system of PDEs. The third array dimension represents the
        % time-step index for a system of time-dependent PDEs, or the
        % eigenvect index for an eigenvalue problem involving a system of
        % PDEs.
        NodalSolution;
        
        % SolutionTimes - Times at which the NodalSolution was computed
        SolutionTimes;
        
        % XGradients - Spatial gradient of solution along x-direction. The
        % shape of the XGradients array is identical to NodalSolution.
        XGradients;
        % YGradients - Spatial gradient of solution along y-direction. The
        % shape of the YGradients array is identical to NodalSolution.
        YGradients;
        % ZGradients - Spatial gradient of solution along z-direction. The
        % shape of the ZGradients array is identical to NodalSolution.
        ZGradients;
        
        
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
        
    end
    
    methods
        function obj = HHXTTransientResults(varargin)
            obj@pde.PDEResults(varargin{:});
            if nargin == 0
                %Will result in a default empty object.
                return
            end
            
            if isempty(obj.NumTimeEig)
                error(message('pde:PDEResults:solNotFromTimeDependent'))
            end
            
            narginchk(3,4);
            pdem = varargin{1};
            u = varargin{2};
            u(isnan(u)) = 0; % maybe add a NaN check here?
            tlist = varargin{3};
            
            if ( ~isnumeric(tlist) || issparse(tlist) || ~isvector(tlist)...
                    || isscalar(tlist) || any(isnan(tlist)) || ...
                    ~any(isfinite(tlist)) )
                error(message('pde:PDEResults:invalidTimeVector'))
            end
            
            t0 = tlist(1);
            tfinal = tlist(end);
            if(t0 == tfinal)
                error(message('pde:PDEResults:tlistEndpointsNotDistinct'));
            end
            
            tdir = sign(tfinal - t0);
            if any( tdir*diff(tlist) <= 0 )
                error(message('pde:PDEResults:tlistUnordered'));
            end
            
            obj.SolutionTimes = tlist;
            
            ureshaped = pde.PDEResults.reshapePDESolution(u, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            
            % For second-order time-dependent PDE(s), nodal derivative is
            % passed as 4th argument.
            if nargin == 4 
                dudt = varargin{4};
                dudtreshaped = pde.PDEResults.reshapePDESolution(dudt, ...
                    obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            else
                dudtreshaped =[];
            end

            if numel(obj.SolutionTimes) ~= obj.NumTimeEig
                error(message('pde:PDEResults:solTimesResultMismatch'))
            end
            obj.NodalSolution = ureshaped;
            obj.NodalTimeDerivative =  dudtreshaped;
            
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
            
            % Calculate gradients, assign nodal gradients as properties and
            % construct interpolants for gradients
            [ux,uy,uz] = nodalGradients(obj.Interpolant);
            obj.XGradients = pde.PDEResults.reshapePDESolution(ux, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            obj.YGradients = pde.PDEResults.reshapePDESolution(uy, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            obj.ZGradients = pde.PDEResults.reshapePDESolution(uz, ...
                obj.IsTimeEig, pdem.PDESystemSize, obj.NumTimeEig);
            
                        
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
            obj.InterpolantXChannelVelocity   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
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
                obj.InterpolantXChannelVelocity   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZChannelVelocity, nstream, obj.NumTimeEig);
                obj.InterpolantZReynolds   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                    obj.ZReynolds, nstream, obj.NumTimeEig);
            end
            
            obj.InterpolantStreamVolHeating   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.StreamVolHeating, nstream, obj.NumTimeEig);
            
            obj.InterpolantStreamVolResistance   = pde.PDEResults.constructInterpolat(pdem.Mesh,...
                obj.StreamVolResistance, nstream, obj.NumTimeEig);
            
            
        end
        
        % Methods declaration
        obj = calcGradients(obj,hhxtm)
        obj = calcVariables(obj,hhxtm)
        
        % Methods declaration
        uintrp = interpolateSolution(obj,varargin)
        [dudxInt, dudyInt, dudzInt] = evaluateGradient(obj,varargin)
        [cdudx, cdudy, cdudz] = evaluateCGradient(obj,varargin)
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
    
    methods (Hidden = true, Static = true, Access = protected)
        function obj = loadobj(obj)
            % function called during loading an object of this type from
            % a MAT-file, interpolant objects need to be constructed
            obj.Interpolant = pde.PDEResults.constructInterpolat(obj.Mesh, ...
                obj.NodalSolution, obj.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdx   = pde.PDEResults.constructInterpolat(obj.Mesh,...
                obj.XGradients, obj.PDESystemSize, obj.NumTimeEig);
            
            obj.InterpolantdUdy   = pde.PDEResults.constructInterpolat(obj.Mesh,...
                obj.YGradients, obj.PDESystemSize, obj.NumTimeEig);
            
            if ~obj.IsTwoD
                obj.InterpolantdUdz   = pde.PDEResults.constructInterpolat(obj.Mesh,...
                    obj.ZGradients, obj.PDESystemSize, obj.NumTimeEig);
            end
        end
        
    end
    
    properties (Hidden = true, SetAccess=protected)
        NodalTimeDerivative
    end
    
    
    
end