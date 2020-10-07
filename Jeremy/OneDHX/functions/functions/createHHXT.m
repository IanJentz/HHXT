function pdem = createHHXT(varargin)
%CREATEPDE - Creates an analysis model 
% PDEM = CREATEPDE(PDESystemSize) creates an empty PDE analysis model that
% can be used to define a system of partial differential equations of size
% PDESystemSize. If a PDESystemSize is not specified, the default size of 1
% is used, giving a scalar system. The output, PDEM, is a PDEModel - a
% container object that provides functions that support the model creation
% through the definition of the geometry, generation of the mesh,
% application of boundary conditions and initial conditions, and solution.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Create a Thermal Analysis Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% PDEM = CREATEPDE('thermal',AnalysisType) creates an empty thermal
% analysis model for the analysis type specified using the AnalysisType
% string. Default AnalysisType is 'steadystate'.
%
%  Thermal AnalysisType         Description    
%-----------------------------------------------------------------
%     'steadystate'        Creates a thermal model for steadystate 
%                          analysis.
%     'transient'          Creates a thermal model for transient 
%                          analysis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Create a Structural Analysis Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% PDEM = CREATEPDE('structural',AnalysisType) creates a structural analysis
% model, for solving small-strain linear elasticity problem. The structural
% model can be used for an analysis of the type specified using
% AnalysisType string from the table below.
%
% Structural AnalysisType         Description    
%-----------------------------------------------------------------
% 'static-solid'          Creates a structural model for static analysis
%                         of a solid (3-D) problem.
% 'static-planestress'    Creates a structural model for static analysis
%                         of a plane-stress problem.
% 'static-planestrain'    Creates a structural model for static analysis 
%                         of a plane-strain problem.
%
% 'modal-solid'           Creates a structural model for modal analysis of 
%                         a solid (3-D) problem.
% 'modal-planestress'     Creates a structural model for modal analysis of 
%                         a plane-stress problem.
% 'modal-planestrain'     Creates a structural model for modal analysis of 
%                         a plane-strain problem.
%
% 'transient-solid'       Creates a structural model for transient dynamics
%                         analysis of a solid (3-D) problem.
% 'transient-planestress' Creates a structural model for transient dynamics
%                         analysis of a plane-stress problem.
% 'transient-planestrain' Creates a structural model for transient dynamics
%                         analysis of a plane-strain problem.
%
% Example 1:  Create a scalar PDEModel, import geometry, assign a
% boundary condition and generate a mesh.
% 
%    numberOfPDEs = 1;
%    pdem = createpde(numberOfPDEs);
%    gm = importGeometry(pdem,'Block.stl')
%    pdegplot(gm, 'FaceLabels','on')
%    % Alternatively, pdegplot(pdem,...) plots the geometry
% 
%    % Assign a boundary condition
%    bc1 = applyBoundaryCondition(pdem,'dirichlet','Face',2, 'u', 0);
% 
%    % Generate and plot the mesh
%    msh = generateMesh(pdem);
%    figure;
%    pdemesh(msh);
%    % Alternatively, pdemesh(pdem,...) plots the mesh
%
% Example 2: Create a thermal model for steady-state analysis.
%
%     thermalmodel = createpde('thermal','steadystate');
%     gm = importGeometry(thermalmodel,'SquareBeam.STL');
%     thermalProperties(thermalmodel,'ThermalConductivity',0.08);
%
% Example 3: Create a structural model for solving a static plane-stress problem.
%
%     structuralmodel = createpde('structural','static-planestress')
%     gm = geometryFromEdges(structuralmodel,@squareg)
%     generateMesh(structuralmodel);
% 
%  See also pde.PDEModel, pde.PDEModel/importGeometry, 
%           pde.PDEModel/applyBoundaryCondition, pde.PDEModel/generateMesh
%           pde.ThermalModel, pde.ThermalModel/thermalProperties
%           pde.StructuralModel, pde.StructuralModel/structuralProperties
%           pde.StructuralModel/structuralBoundaryLoad, pde.StructuralModel/structuralBC

%  Copyright 2014-2017 The MathWorks, Inc.

    if nargin > 0
        [varargin{:}] = convertStringsToChars(varargin{:});
    end

    narginchk(0,2);
    
    if (nargin == 0)
        pdem = pde.PDEModel(1);
%         pdem = hhxt.HHXTModel(1);
        addprop(pdem,'MaterialProperties');
        addprop(pdem,'MaterialMap');
        addprop(pdem,'AnalysisType');
        return
    end
    
    if (nargin == 1) && isnumeric(varargin{1})
%         validateattributes(varargin{1},{'numeric'},{'real', 'finite', 'integer', 'nonsparse', 'nonnan'});
%         pdem = pde.PDEModel(double(varargin{1}));
        pdem = hhxt.HHXTModel(double(varargin{1}));
%         addprop(pdem,'MaterialProperties');
%         addprop(pdem,'MaterialMap');
%         addprop(pdem,'AnalysisType');
%         addprop(pdem,'MatrixOptions');
%         addprop(pdem,'WorkingDirectory');
%         addprop(pdem,'PermeabilityFloor');
%         addprop(pdem,'AdvectTermDemotion');   
%         addprop(pdem,'AddtEquationCoefficients');
        return
    end
   
    
    if (nargin > 1) && isnumeric(varargin{1})
        error(message('pde:createpde:invalidFirstArg'));
    end
    
     if isempty(varargin{1}) % String special cases: '',"",0x1,0x0, as convertStringsToChars returns an empty cell for 0x1,0x0 strings.
         error(message('pde:createpde:unsupportedVertical',''));
     end
     
    %Only valid input argument from this point on must be a char/string.
    validateattributes(varargin{1}, {'char','string'},{'scalartext'});
   


end

