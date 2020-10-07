classdef (Sealed) HHXTMaterialAssignment < handle & matlab.mixin.internal.CompactDisplay & matlab.mixin.Copyable
% HHXTStreamMaterialAssignment - Specify stream material properties over domain or subdomain
%  For a HHXT model, the toolbox solves equations of the following form:
%
%  (MassDensity*SpecificHeat)*dT/dt - div(ThermalConductivity*grad(T)) = HeatSource
%
%     This method creates an object representing the fluid properties:
%     MassDensity, SpecificHeat, ThermalConductivity, and Viscosity; 
%     and Channel properties: HydraulicDiameter, fDarcy, and jColburn; of
%     a stream in a domain or subdomain and appends the object to the 
%     MaterialProperties property.
%
%     Instances of this class can only be created by calling the
%     HHXTProperties method of the HHXTModel class.
%
% See also pde.ThermalModel, pde.ThermalModel/thermalProperties

% Copyright 2018 Ian Jentz
    
properties (SetAccess = private)
    % RegionType - Type of geometric region the coefficients are assigned to.
    %    A string specifying the type of geometric domain the coefficients
    %    are assigned to. This string has two possible values: 'cell'
    %    or 'face'.
    RegionType;
end
    
properties
    % RegionID - ID of the geometric regions the coefficients are assigned to.
    %    For 2-D each RegionID satisfies 0 < RegionID(j) < NumFaces in the
    %    geometry. For 3-D each RegionID satisfies 0 < RegionID(j) < NumCells
    %    in the geometry.
    RegionID;
    
    % StreamID - ID of the stream region the coefficients are assigned to.
    %    For the HHXT model StreamID satisfies 0 < StreamID(j) < NumStreams
    %    where StreamID = 0 corresponds to the core solid material and
    %    StreamID(j) > 1 corresponds to fluid streams
    StreamID;
    
    % PhaseFraction - Volume fraction that the fluid occupies. 
    %    For the solid features this corresponds to the volume fraction
    %    of material occupied in the core.  For stream features this
    %    corresonds to the porosity of the core.
    PhaseFraction;

    % ThermalConductivity of the material
    %   A numerical scalar or a function handle that
    %   computes thermal conductivity for a give spatial location at a
    %   given temperature and pressure.
    ThermalConductivity;

    % MassDensity of the material
    %   A numerical scalar or a function handle that
    %   computes mass density for a give spatial location at a
    %   given temperature and pressure.
    MassDensity;

    % SpecificHeat of the material
    %   A numerical scalar or a function handle that
    %   computes specific heat capacity for a give spatial location at a
    %   given temperature and pressure.
    SpecificHeat;
    
    % Enthalpy of the fluid
    %   A numerical scalar or a function handle that
    %   computes the enthalpy for a give spatial location at a
    %   given temperature and pressure.
    Enthalpy;
    
    % Viscosity of the fluid
    %   A numerical scalar or a function handle that
    %   computes viscosity for a given spatial location at a given
    %   temperature and pressure.
    Viscosity;
    
    % HydraulicDiamter
    %   A numerical scalar or a function handle that
    %   computes hydraulic diameter for a given spatial location.
    HydraulicDiameter;
    
    % fDarcy
    %   A numerical scalar or matrix, or a function handle that 
    %   computes the Darcy friction factor for a given spatial location
    %   at a given Reynolds number.  Note Reynolds number will have to be
    %   defined in the function based on velocities and viscosity.
    fDarcy;
    
    % PermeabilityTransform
    %   A transformation matrix that operates on the permiability matrix.
    %   This can be used to have spatially variying permeability within a
    %   region.
    PermeabilityTransform;

    % Permeability
    %   A fixed Darcy Porous media permiability.
    Permeability;
    
    % NodalPermeability
    %   Permiability values, k_xx, k_yy, k_zz defined at nodes
    %   and stored as pde results for later interpolation
    NodalPermeability
    
    % NodalDarcyVelocity
    %   Darcy Velocity values, uD_z, uD_y, and uD_z defined at nodes
    %   and stored as pde results for later interpolation
    NodalDarcyVelocity
    
    % NodalDarcyVelocity
    %   Temperature avdection values, uD_z*dT/dx, uD_y*dT/dy, and uD_z*dT/dz defined at nodes
    %   and stored as pde results for later interpolation
    NodalTempAdvection
    
    % jColburn
    %   A numerical scalar or a function handle that 
    %   computes the Colburn heattransfer coefficient for a given spatial 
    %   location at a given Reynolds number.
    jColburn;
    
    % Nusselt
    %   A numerical scalar or a function handle that 
    %   computes the Colburn heattransfer coefficient for a given spatial 
    %   location at a given Reynolds number.
    Nusselt;
    
    % WallFraction
    %    Fraction or the core solid that is assigned to the wall. 
    %    A portion of the core solid mass is assigned to the fluid wall.
    %    This specified as the fraction of core mass to be assigned to the
    %    wall.
    WallFraction;
    
    % WallCoreResistance
    %    Characteristic Volumetric Thermal resistance [W-m3/K].
    Resistance;
    
end
    
methods %ctor
    function self = HHXTMaterialAssignment(tmar,varargin)
        self.RecordOwner = tmar;
        parser = inputParser;
        parser.addParameter('face', []);
        parser.addParameter('cell', []);
        parser.addParameter('stream', []);
        parser.addParameter('PhaseFraction',[]);
        parser.addParameter('ThermalConductivity', []);
        parser.addParameter('MassDensity', []);
        parser.addParameter('SpecificHeat', []);
        parser.addParameter('Enthalpy', []);
        parser.addParameter('Viscosity', []);
        parser.addParameter('HydraulicDiameter', []);
        parser.addParameter('fDarcy', []);
        parser.addParameter('PermeabilityTransform', []);
        parser.addParameter('jColburn', []);
        parser.addParameter('Nusselt', []);
        parser.addParameter('WallFraction', []);
        parser.addParameter('Resistance', []);
        parser.addParameter('NodalPermeability',[]);
        parser.addParameter('NodalDarcyVelocity',[]);
        parser.addParameter('NodalTempAdvection',[]);
        parser.addParameter('Permeability', []);
        parser.addParameter('Material',[]); % previously defined HHXTMaterialAssignment
        parser.parse(varargin{:});
        
        if ~isempty(parser.Results.Material)
            self = copy(parser.Results.Material);  % to use make sure to have matlab.mixin.Copyable subclassed
        end

        numdims = 2;
        if ~isempty(parser.Results.face)
            self.RegionType = 'face';
            self.RegionID = parser.Results.face;
        elseif ~isempty(parser.Results.cell)
            self.RegionType = 'cell';
            self.RegionID = parser.Results.cell;
            numdims = 3;
        end
        systemsize = 1;
        
        if ~isempty(parser.Results.stream)
            self.StreamID = parser.Results.stream;        
        elseif isempty(parser.Results.Material)
            self.StreamID = 0;
        end
        if ~isempty(parser.Results.PhaseFraction)
%             validateattributes(parser.Results.PhaseFraction,{'double'},{'scalar','>=',0,'<=',1});
            self.PhaseFraction = parser.Results.PhaseFraction;
        end
        
        if ~isempty(parser.Results.ThermalConductivity)
            self.ThermalConductivity = parser.Results.ThermalConductivity;
        end
        if ~isempty(parser.Results.MassDensity)
            self.MassDensity = parser.Results.MassDensity;
        end

        if ~isempty(parser.Results.SpecificHeat)
            self.SpecificHeat = parser.Results.SpecificHeat;
        end
        
        if ~isempty(parser.Results.Enthalpy)
            self.Enthalpy = parser.Results.Enthalpy;
        end
        
        if ~isempty(parser.Results.Viscosity)
            self.Viscosity = parser.Results.Viscosity;
        end
        
        if ~isempty(parser.Results.HydraulicDiameter)
            self.HydraulicDiameter = parser.Results.HydraulicDiameter;
        end
        
        if ~isempty(parser.Results.fDarcy)
            self.fDarcy = parser.Results.fDarcy;
        end
        
        if ~isempty(parser.Results.PermeabilityTransform)
            self.PermeabilityTransform = parser.Results.PermeabilityTransform;
        end
        
        if ~isempty(parser.Results.jColburn)
            self.jColburn = parser.Results.jColburn;
        end
        
        if ~isempty(parser.Results.Nusselt)
            self.Nusselt = parser.Results.Nusselt;
        end
        
        if ~isempty(parser.Results.WallFraction)
            self.WallFraction = parser.Results.WallFraction;
        end
        
        if ~isempty(parser.Results.Resistance)
            self.Resistance = parser.Results.Resistance;
        end

        if ~isempty(parser.Results.Permeability)
            self.Permeability = parser.Results.Permeability;
        end

        if ~isempty(parser.Results.NodalPermeability)
            self.NodalPermeability = parser.Results.NodalPermeability;
        end
        
        if ~isempty(parser.Results.NodalDarcyVelocity)
            self.NodalDarcyVelocity = parser.Results.DarcyVelocity;
        end
        
        if ~isempty(parser.Results.NodalTempAdvection)
            self.NodalTempAdvection = parser.Results.NodalTempAdvection;
        end


%         These checks can be put in later
%         self.checkAllMatrixMtlSizes(systemsize, numdims);
%         self.checkFcnHdlArgCounts(systemsize, numdims);
    end
end


methods % Setter methods

    function set.RegionType(self, rtype)
        self.RegionType = rtype;
    end

    function set.RegionID(self, rids)
        self.ValidateRegionID(rids);
        self.RegionID = rids;
    end
    
    function set.StreamID(self, sids)
        self.ValidateStreamID(sids);
        self.StreamID = sids;
    end
    
    function set.PhaseFraction(self, coef)
        self.ValidatePhaseFraction(coef,'PhaseFraction');
        self.PhaseFraction = coef;
    end

    function set.MassDensity(self, coef)
        self.CoefPrecheck(coef,'MassDensity');
        self.MassDensity = coef;
    end

    function set.SpecificHeat(self, coef)
        self.CoefPrecheck(coef,'SpecificHeat');
        self.SpecificHeat = coef;
    end
    
    function set.Enthalpy(self, coef)
        self.CoefPrecheck(coef,'Enthalpy');
        self.Enthalpy= coef;
    end


    function set.ThermalConductivity(self, coef)
        self.CoefPrecheck(coef,'ThermalConductivity');
        self.ThermalConductivity = coef;
    end

    function set.Viscosity(self, coef)
        self.CoefPrecheck(coef,'Viscosity');
        self.Viscosity = coef;
    end
    
    function set.HydraulicDiameter(self, coef)
        self.CoefPrecheck(coef,'HydraulicDiameter');
        self.HydraulicDiameter = coef;
    end
    
    function set.fDarcy(self, coef)
        self.CoefPrecheck(coef,'fDarcy');
        self.fDarcy = coef;
    end
    
    function set.PermeabilityTransform(self, coef)
        self.ValidatePermeabilityTransform(coef,'PermeabilityTransform');
        self.PermeabilityTransform = coef;
    end
    
    function set.jColburn(self, coef)
        self.CoefPrecheck(coef,'jColburn');
        self.jColburn = coef;
    end
    
    function set.Nusselt(self, coef)
        self.CoefPrecheck(coef,'Nusselt');
        self.Nusselt = coef;
    end
    
    function set.WallFraction(self, coef)
        self.CoefPrecheck(coef,'WallFraction');
        self.WallFraction = coef;
    end
    
    function set.Resistance(self, coef)
        self.CoefPrecheck(coef,'Resistance');
        self.Resistance = coef;
    end

    function set.Permeability(self, coef)
        self.ValidatePermeability(coef,'Permeability');
        self.Permeability = coef;
    end
    
    function set.NodalPermeability(self, coef)
        self.ValidateNodalResults(coef,'NodalPermeability');
        self.NodalPermeability = coef;
    end
    
    function set.NodalDarcyVelocity(self, coef)
        self.ValidateNodalResults(coef,'NodalDarcyVelocity');
        self.NodalDarcyVelocity = coef;
    end
    
    function set.NodalTempAdvection(self, coef)
        self.ValidateNodalResults(coef,'NodalTempAdvection');
        self.NodalTempAdvection = coef;
    end



end


methods(Static, Access = private)
    function tf = coefDefined(coef)
        tf = (isnumeric(coef) && ~(isscalar(coef) && coef == 0) || isa(coef,'function_handle'));
    end

    
    function ok=ValidateRegionID(rgnid)
        % Must be real(non-complex), full, natural number.
        if ~isreal(rgnid) || ~all(rgnid(:) > 0) || issparse(rgnid) || any(mod(rgnid(:),1)~=0)
            error(message('invalid Region ID'));
        end
        ok = true;
    end
    
    function ok=ValidateStreamID(strid)
        % Must be real(non-complex), full, natural number.
        if ~isreal(strid) || ~all(strid(:) >= 0) || issparse(strid) || any(mod(strid(:),1)~=0)
            error(message('invalid Stream ID'));
        end
        ok = true;
    end
    
    function ok=ValidatePhaseFraction(coef,coefName)
        % Must be real(non-complex), full, natural number.
%         if ~isreal(phase) || ~all(phase(:) >= 0) || issparse(phase) || ~all(phase(:) <= 1)
%             error(message('invalid Phase Fraction'));
%         end
        validateattributes(coef,{'double','function_handle','cell'},{},'',coefName)
        ok = true;
    end
    
    function ok=ValidatePermeability(coef,coefName)
        % Must be a double, or cell
        validateattributes(coef,{'double','cell'},{},'',coefName)
        
        ok = true;
    end
    
    function ok=ValidatePermeabilityTransform(coef,coefName)
        % Must be a double, function, or cell
        validateattributes(coef,{'double','function_handle','cell'},{},'',coefName)
        
        ok = true;
    end
    
    function ok=ValidateNodalResults(coef,coefName)
        % must be a steady state or transient pde result
        validateattributes(coef,{'pde.StationaryResults','pde.TransientResults','double'},{},'',coefName)
        
        ok = true;
    end
    
    function ok=CoefPrecheck(coef,coefName)
        
      validateattributes(coef,{'double','function_handle'},{},'',coefName)
        if isa(coef,'function_handle')
            validateattributes(coef,{'function_handle'},{'scalar'},'',coefName);
        else
            validateattributes(coef,{'double'},{'nonnan','nonempty','nonsparse','real','finite','nonnegative'},coefName);
            if all(coef==0) % Matrix form of coef can have some zeros, but not all.
                error(message('Invalid zero valued material property'));
            end
        end

        ok = true;
    end
    
end


% methods(Hidden=true)
%     function delete(self)
%         if ~isempty(self.RecordOwner) && isvalid(self.RecordOwner)
%             self.RecordOwner.delistMaterialAssignments(self);
%         end
%     end
% end



properties (Hidden = true, Access='private')
    RecordOwner;
end




end

