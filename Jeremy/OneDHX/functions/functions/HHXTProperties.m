function mtl = HHXTProperties(self,varargin)
%HHXTProperties - Assign properties of the HHXT multistream to a HHXT model
%   This acts as a construction function for a HHXTMaterialAssignment.
%   Parsing of inputs is performed and geometric specification is checked
%   against dimensionality of the model.
%
%   mtl = HHXTProperties(HHXTmodel,PROPERTY,VALUE) assigns material
%   properties to a HHXT model.  VALUE can be specified as
%       (1) A numerical scalar for constant material properties
%       (2) MATLAB function format for a non-constant or nonlinear material
%
%   PROPERTY   VALUE
%--------------------------------------------------------------------------
%   'PhaseFraction'
%   'ThermalConductivity'
%   'MassDensity'
%   'SpecificHeat'
%   'HydraulicDiameter'
%   'Viscosity'
%   'jColburn'
%   'fDarcy'
%   'WallFraction'
%   'Resistance'
%
%  For steady-state analysis the only required material properties are:
%      ThermalConductivity, in the solid core (stream = 0)
%      PhaseFraction, MassDensity, SpecificHeat, Viscosity,
%    HydraulicDiameter, jColburn, fDarcy, and WallCoreResistance, in the
%    fluid streams (stream > 0)
%  For transient analysis the following material properties are required:
%      PhaseFraction, ThermalConductivity, MassDensity, and SpecificHeat,
%    in the solid core (stream = 0)
%      PhaseFraction, MassDensity, SpecificHeat, Viscosity,
%    HydraulicDiameter, jColburn, fDarcy, WallFraction, and WallCoreResistance,
%    in the fluid streams (stream > 0)
%
%  mtl = thermalProperties(__,'Face',FACEID) assigns properties for
%  the specifed face of a 2-D geometry. FACEID is a positive integer in the
%  range of 1 to thermalmodel.Geometry.NumFaces.
%
%  mtl = thermalProperties(__,'Cell',CELLID) assigns properties for
%  the specifed cell of a 3-D geometry. CELLID is a positive integer in the
%  range of 1 to thermalmodel.Geometry.NumCells.
%
%  mtl = thermalProperties(__,'Stream',STREAMID) assigns properties for
%  the specifed stream. CELLID is a positive integer in the
%  range of 0 to number of streams.

narginchk(3,27);
nargoutchk(0,1);

if isempty(self.Geometry)
    error('HHXT Model Geometry not specified');
end

% determine the number of streams
N = self.PDESystemSize;
numStreams = 0;
if ~mod((N-1),2) % multiple streams no wall participation
    numStreams = (N-1)/2;
elseif ~mod((N-1),3) % multiple streams with wall participation
    numStreams = (N-1)/3;
end

parser = inputParser;
parser.addParameter('face', [], @pde.EquationModel.isValidEntityID);
parser.addParameter('cell', [], @pde.EquationModel.isValidEntityID);
parser.addParameter('stream', [], @isValidStreamAttribute);
parser.addParameter('MassDensity', [],@isValidScalarAttribute);
parser.addParameter('SpecificHeat', [],@isValidScalarAttribute);
parser.addParameter('Enthalpy', [],@isValidScalarAttribute);
parser.addParameter('ThermalConductivity', [],@isValidMatrixAttribute);
parser.addParameter('PhaseFraction', [],@isValidScalarAttribute);
parser.addParameter('HydraulicDiameter', [],@isValidScalarAttribute);
parser.addParameter('Viscosity', [],@isValidScalarAttribute);
parser.addParameter('jColburn', [],@isValidScalarAttribute);
parser.addParameter('Nusselt', [],@isValidScalarAttribute);
parser.addParameter('fDarcy', [],@isValidMatrixAttribute);
parser.addParameter('PermeabilityTransform', [],@isValidMatrixAttribute);
parser.addParameter('WallFraction', [],@isValidScalarAttribute);
parser.addParameter('Resistance', [],@isValidScalarAttribute);
parser.addParameter('Permeability', [],@isValidMatrixAttribute);
parser.addParameter('NodalPermeability', [],@isValidPDEResults);
parser.addParameter('NodalDarcyVelocity', [],@isValidPDEResults);
parser.addParameter('NodalTempAdvection', [],@isValidPDEResults);
parser.addParameter('Material',[],@isValidMaterialAttribute);

parser.parse(varargin{:});
argsToPass = varargin;

%If a face (for 2-D) or a cell (for 3-D) is not explicitly specified then
%use the same material for all faces and edges.
if self.IsTwoD % problem is 2D
    if ~isempty(parser.Results.cell)
        error('cell properties cannot be specified in a 2D analysis');
    end
    if isempty(parser.Results.face)
        argsToPass = [{'face',1:self.Geometry.NumFaces}, varargin];
    else
        self.isValidEntityID(parser.Results.face);
        if any(parser.Results.face > self.Geometry.NumFaces)
            error(message('invalid face index'));
        end
    end
else % problem is 3D
    if ~isempty(parser.Results.face)
        error(message('face properties cannot be specified in a 3D analysis'));
    end    
    if isempty(parser.Results.cell)
        argsToPass = [{'cell',1:self.Geometry.NumCells}, varargin];
    else
        self.isValidEntityID(parser.Results.cell);
        if any(parser.Results.cell > self.Geometry.NumCells)
            error(message('invalid cell index'));
        end
    end
    
end
% check stream input
stream = 0;
if (isempty(parser.Results.stream) && isempty(parser.Results.Material))
    argsToPass = [{'stream',0},varargin]; 
elseif isempty(parser.Results.Material)
    isValidStreamAttribute(parser.Results.stream);
    if any(parser.Results.stream > numStreams)
        error('invalid stream index')
    end
    stream = parser.Results.stream;
end

% if parser.Results.stream == 0 % adding solid thermal properties
% else %adding fluid stream properties
% end

if isempty(self.MaterialProperties)
    mtlContainer = pde.MaterialAssignmentRecords(self);
%         if stream == 0
%             mtl = HHXTSolidMaterialAssignment(mtlContainer,argsToPass{:});
%         else
%             mtl = HHXTStreamMaterialAssignment(mtlContainer,argsToPass{:});
%         end
    mtl = hhxt.HHXTMaterialAssignment(mtlContainer,argsToPass{:});
    mtlContainer.MaterialAssignments = mtl;
    self.MaterialProperties = mtlContainer;
else
%         if stream == 0
%             mtl = HHXTSolidMaterialAssignment(self.MaterialProperties,argsToPass{:});
%         else
%             mtl = HHXTStreamMaterialAssignment(self.MaterialProperties,argsToPass{:});
%         end
    mtl = hhxt.HHXTMaterialAssignment(self.MaterialProperties,argsToPass{:});
    self.MaterialProperties.MaterialAssignments(end+1) = mtl;
end

%     mtlContainer = pde.MaterialAssignmentRecords(self);
%     if stream == 0
%         mtl = HHXTSolidMaterialAssignment(mtlContainer,argsToPass{:});
%     else
%         mtl = HHXTStreamMaterialAssignment(mtlContainer,argsToPass{:});
%     end


    function valid = isValidStreamAttribute(x)
    % determines if scalar inputs are valid
    % can be: doubles | function handle
        validateattributes(x,{'double'},{'scalar','integer','>=',0});
        valid = true;
    end

    function valid = isValidMaterialAttribute(x)
        validateattributes(x,{'hhxt.HHXTMaterialAssignment','HHXTMaterialAssignment'},{})
        valid = true;
    end

    function valid = isValidScalarAttribute(x)
    % determines if scalar inputs are valid
    % can be: doubles | function handle
        validateattributes(x,{'double','function_handle'},{'scalar'});
        valid = true;
    end

     function valid = isValidMatrixAttribute(x)
    % determines if matrix inputs are valid
    % can be: doubles | function handle
        validateattributes(x,{'double','function_handle'},{});
        valid = true;
     end

    function valid = isValidPDEResults(x)
        validateattributes(x,{'pde.StationaryResults','pde.TransientResults'},{});
        valid = ture;
    end



end

