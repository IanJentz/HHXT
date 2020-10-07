function bc = applyStreamBoundaryCondition(self,varargin)
% applyStreamBoundaryCondition - Apply a boundary condition (BC) to the
% HHXT model geometry.  Creates a BC object specific to the HHXT model that
% applies a BC to the geometry stored in the Geometry property.  The BC
% object is appended to the StreamBoundaryConditionsRecords and a handle to
% the BC is returned.
%
% The HHXT toolbox also supports BCs specified using the PDE toolbox's
% applyBounaryCondition().
%
% BC = applyStreamBoundaryCondition(HHXTM,...) sets the boundary conditions
% of model HHXTM.  Inputs are specified at name-value pairs following the
% model input.
%
% BC = applyStreamBoundaryCondition(...,REGIONTYPE, REGIONID,...) sets the
% boundary conditions on a region of the domain defined by REGIONTYPE and
% REGIONID. Where name/value pari REGIONTYPE is 'Face','Edge', 'Vertex' and
% REGIONID is an integer specifying the ID of the geometric entity.
%
% BC = applyStreamBoundaryCondition(...,'Stream',STREAMID,...) sets the
% boundary condition to apply to the stream defined by STREAMID.  Where
% STREAMID is an integer specifying the ID of the HHXT model stream, 
% 0 for the solid of the HHXT model and 1,2,...(N-1)/2 for the various streams.
%
% BC = applyStreamBoundaryCondition(...,'Direction',FLOWDIR,...) sets the
% boundary condition to be an inlet or an outlet as defined by FLOWDIR.
% Where FLOWDIR is a string specifying either 'inlet' or 'outlet'.  The
% flow direction need not be specified for BCs applying to the HHXT
% solid (STREAMID = 0).
%
% BC = applyStreamBoundaryCondition(...,'name','value') sets the boundary
% conditions using the following name-value pairs.
%
% For Dirichlet boundary conditions, specify 'Temperature' and/or 'Pressure'
%       'Temperature'       temperature in C of the BC.  This is the solid
%                           temperature if STREAMID = 0, otherwise it is
%                           the fluid temperature.  Note specifying the
%                           fluid temperature of an outlet is not
%                           recommended as the fluid will be flowing out of
%                           the domain, leaving only fluid conductivity,
%                           and no advection, to fix the BC.
%       'Pressure'          pressure in Pa of the BC.  This does not apply
%                           to the solid when STREAMID = 0, otherwise it is
%                           the fluids stagnation pressure.
%
% For Neuman boundary conditions, specify
% 'HeatFlux','HeatFlow','MassFlux', or 'MassFlow' conditions using the
% following name-value pairs.
%       'HeatFlux'          heat flux in W/m2 through the BC.
%
% See also pde.BoundaryCondition, pde.PDEModel/geometryFromEdges, DECSG

%   Detailed explanation goes here

    narginchk(3,14); 
    nargoutchk(0,1);
    
    argsToPass = [varargin,{'SystemSize', self.PDESystemSize}];
    
    nc = numel(varargin{1});
    if ~(strncmpi(varargin{1},'stream',nc) || strncmpi(varargin{1},'face',nc) || strncmpi(varargin{1},'edge',nc) || strncmpi(varargin{1},'cell',nc)...
            || strncmpi(varargin{1},'Temperature',nc) || strncmpi(varargin{1},'HeatFlux',nc) || strncmpi(varargin{1},'HeatFlow',nc) || strncmpi(varargin{1},'VolumetricGeneration',nc)...
            || strncmpi(varargin{1},'Pressure',nc) || strncmpi(varargin{1},'MassFlux',nc) || strncmpi(varargin{1},'MassFlow',nc))
        % not a valid first argument
        error(message('pde:pdeBoundaryConditions:invalidBCType'));
    end
    
    argsToTest = argsToPass(1:end);
    
    if ~isempty(self.Geometry) && numel(varargin) > 1  
        parser = inputParser;                   
        parser.KeepUnmatched=true;  
        parser.addParameter('Face', [], @isreal);
        parser.addParameter('Edge', [], @isreal);
        parser.addParameter('Stream', []);
        parser.addParameter('Direction', []);
        parser.addParameter('Temperature', []);
        parser.addParameter('HeatFlux', []);
        parser.addParameter('HeatFlow', []);
        parser.addParameter('Pressure', []);
        parser.addParameter('MassFlux', []);
        parser.addParameter('MassFlow', []);
        parser.parse(argsToTest{:}); 
        
        if self.IsTwoD && ~isempty(parser.Results.Face)
            error('a 2D problem can only have a boundary condition on an Edge')
        end
        if ~self.IsTwoD && ~isempty(parser.Results.Edge)
            error('a 3D problem can only have a boundary condition on an Face')
        end
        
        if ~isempty(parser.Results.Face) 
            validateattributes(parser.Results.Face,{'numeric'},{'integer','positive', 'nonzero','real', 'nonsparse'});           
            if any(parser.Results.Face > self.Geometry.NumFaces)
              error(message('pde:pdeModel:invalidFaceIndex'));
            end
        elseif ~isempty(parser.Results.Edge)           
            validateattributes(parser.Results.Edge,{'numeric'},{'integer','positive','nonzero','real','nonsparse'});           
            if any(parser.Results.Edge > self.Geometry.NumEdges)
              error(message('pde:pdeModel:invalidEdgeIndex'));
            end
        end
    
    
        if ( ~isempty(parser.Results.Face) + ~isempty(parser.Results.Edge) ) ~= 1
            error('a Face, or Edge region must be specified')
        end
        
        heatonly = false;
        if ~isempty(parser.Results.Stream)
            if parser.Results.Stream == 0
                heatonly = true;
            end
            n = self.PDESystemSize;
            nstr = (n-1)/2;
            if parser.Results.Stream > nstr
                error(['unrecognized stream, model only contains ',num2str(nstr),' streams']);
            end
        else
            heatonly = true;
        end
        if heatonly && any([~isempty(parser.Results.Pressure),~isempty(parser.Results.MassFlux),~isempty(parser.Results.MassFlow)])
            error('boundary conditions for stream 0 (solid) can only be Temperature, HeatFlux, or HeatFlow')
        end
        if ~heatonly && isempty(parser.Results.Direction)
            error('stream flow direction of inlet or outlet must be specified')
        end
      

%         bc = hhxt.StreamBoundaryCondition(self,varargin{:});
        
        if isempty(self.StreamBoundaryConditions)
            BCcont = hhxt.StreamBoundaryConditionRecords(self);             
            bc = hhxt.StreamBoundaryCondition(BCcont,argsToPass{:});      
            BCcont.BoundaryConditionAssignments = bc; 
            self.StreamBoundaryConditions = BCcont;
        else
            bc = hhxt.StreamBoundaryCondition(self.StreamBoundaryConditions,argsToPass{:});
            self.StreamBoundaryConditions.BoundaryConditionAssignments(end+1) = bc;         
        end  
    
    end

end

