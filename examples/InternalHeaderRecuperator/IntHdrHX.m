%% Define how this script runs

breakpoints = false; %will pause execution at points, press F5 to continue
verbose = true; % will write progress and output to command line
plotting = true; % choose to plot as the script progresses

% choose how to model the deviation of the channels within the transition,
% entrance, and exit regions.
% 'simple' uses a 29 degree rotation of the channels throughout
% 'complex' uses a more realistic mulit-rotational position dependent path
% see InternalHeader_MicrochannelTransforms.pdf for illustration
GammaMode = 'complex';
<<<<<<< Updated upstream

flowBC = 'mdot'; %'DP', or 'mdot' for Pressure drop, or mass flow Boundary Conditions
=======
>>>>>>> Stashed changes

%file handeling
basename = 'IntHdrHX_Results';
filename = [basename,'_solu']; % steady state results filename
file_init = [basename,'_init']; % intial condition results filename

%% HX geometry and condition variables
% These variables define the geometry of the rectangular model, the
% conditions imposed at the boundaries, the microchannel properties, and
% the friction factor and Nusselt number of the microchannels.

%geometry
H_max = 0.005;  %mesh size, maximum size of an element in m^2
H = 77e-3;    %height in m (although the model is 2D it has a finite thickness)

%cold side conditions
m_dot_C = 0.3;      %mass flow in kg/s
T_C_in = 50;        %temperature in C
T_C_out = 592;      %temperature in C
P_C_in = 27.5e6;    %pressure in Pa
DP_C = 0.31e6;       %pressure drop in Pa
P_C_out = P_C_in - DP_C; 

%hot side conditions
m_dot_H = 0.3;      %mass flow in kg/s
T_H_in = 650;       %temperature in C
T_H_out = 85;       %temperature in C
P_H_in = P_C_out;    %pressure in Pa
DP_H = 0.32e6;       %pressure drop in Pa
P_H_out = P_H_in - DP_H;

% we define the micro-channel geometry.  Our HX consists of 27 hot plates
% and 28 cold plates, each containing the same Zig-zag micro-channel
% note the phase fraction varaibels are global so that they can be used
% within functional phase fraction definitions, e.g. @phi_SHC_trans
D_h = 0.9523e-3; % the hydraulic diameter of the micro-channel in m
global phi_C phi_H
phi_core = 0.2041; % volume fraction of channels in the core
phi_C = ((28/55)*phi_core); % volume fraction of cold channels within the core (28 of 55 plates)
phi_H = ((27/55)*phi_core); % volume fraction of hot channels within the core (27 of 55 plates)

%dimensionless channel performance parameters
f_D = 0.181; %Darcy friction factor four times the fanning factor
fmatr = [f_D;100*f_D];
Nu = 59.4;

%% Problem Initialization
% Construct the HHXT model object using the |createHHXT()| function.  Set the
% analyis to be a steady state problem using the AnalysisType property.
% Set the solution options using the MatrixOptions property.  Set the
% solver options using the SolverOptions property

% create an HHXT model consiting of 2 fluid streams (hot and cold)
% this will require a N=5 degree of freedom problem
% the d.o.f's are solid Temp, cold side Press, cold side Temp, hot side
% Press and hot side Temp.
model = createHHXT(5); %constructor for the HHXT environment, N=5
model.WorkingDirectory = cd;
% set the analysis to be steady state
model.AnalysisType = 'steadystate';

% the problem will be solved in this way...
model.MatrixOptions = {...
    'Permeability','calculated',... % with permeability calculated from a Darcy friction factor 
    'ResistanceMode','Nusselt', ... % with Thermal resistance to heat transfer defined by a Nusselt number
    'Conductivity','aniso',...      % with anisotropic conductivity in the solid body
    'FluidConduction','on','Advection','on',... % with both conductive and advective heat transport within the fluid streams
    'PermeabilityTransform','on'};  % we turn on permeability transform as 
                                    % we are defining rotational transformations 
                                    % of the microchannel direction, Gamma

% the solver options are...
model.SolverOptions.MaxIterations = 20;
model.SolverOptions.ResidualTolerance = 1e-6;
model.SolverOptions.MinStep = 1;  % smaller steps (<1) can be taken but are not necessary
model.SolverOptions.AbsoluteTolerance = 1e-2;
if verbose == true
    model.SolverOptions.ReportStatistics = 'on'; % will write iterations and convergence data to the command line
end

%% Geometry
% Define the geometry of the HX using the |Geom_IntHdrHX()| function.
% This will result in a geometry with 10 Faces/Regions and 54 edges
%
% The internal header heat exchanger has fluid inlets and outlets that are
% contained within the recuperator body.  Instead of having fluids exit at
% the sides or end of the heat exchanger the hot and cold sides are
% seperated and routed to different circulare intlet and outlet plenums.
%
% This is accomplished using the following features:
% - core region: face 3, this is the main zig-zag microchannel core
% - jump region: faces 1 and 10, here a pair of zig-zag channels is
% combined into a straight channel, halving the total number of channels
% - transition region: faces 4 and 6, here straight channels are diverged
% to seperate hot and cold sides of the heat exchanger
% - entrace/exit regions: faces 2, 5, 8, and 9, contain only hot or cold
% straight channels that are directed to inlet and outlet boundaries
% - solid wall: face 7, completely solid material that forms the side wall
% and part of the internal header walls
%

G = Geom_IntHdrHX(model);
model.setThickness(H);  % set the 2D thickness (if not set the default will be 1 m thick!)

if plotting == true % plot the geometry in a figure
figure()
pdegplot(model,'CellLabels','on','EdgeLabels','on','FaceLabels','on','FaceAlpha',0.5);
end

if breakpoints == true
disp('....press F5 to continue')

keyboard %programatically inserts a breakpoint
end

%% Mesh the problem
% Mesh the geometry using the |generateMesh()| function.  The HHXT model only
% works with linear elements only at this point.
%
% Refine the mesh once within the internal header regions.
% Use the |refineHHXTmesh()| function.

% generate a mesh of linear tets
if verbose == true
disp('    ...generating mesh')
end
generateMesh(model,'Hmax',H_max,'GeometricOrder','linear');

% calling refineHHXTmesh with a single integer or row of integers
% will refine the mesh within the region/face that is specified.
model = refineHHXTmesh(model,[1,2,4,5,6,8,9,10]);    % refine triangles within Faces 1,2 ,3 ,4 ,6, 8, 9, & 10
% model = refineHHXTmesh(model,[1,2,5,8,9,10]);    % refine triangles within Faces 1,2 ,3, 8, 9, & 10

if plotting == true % plot the geometry in a figure
figure()
pdeplot(model,'FaceAlpha',0.5);  %plot the mesh using the pde toolbox's pdeplot function
end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint
end

%% Material Region Definitions
% Define the behavior of the HHXT model using the |HHXTProperties()|
% function.  This sets the homogenized properties of the HX for each region
% of the geometry.  We have to define conditions in 10 seperate region,
% Faces 1-10.
%
% The internal header recuperator presents two unique challenges:
% 1) The number and distribution of channels change within the jump,
% transition, and entrance and exit regions (faces 1, 2, 4, 5, 6, 8, 9, & 10) 
% The volume fraction, phi, of the fluid streams and solid body
% change within these regions.  This is illustrated in InternalHeader_VolFractionChanges.pdf
% Functional definitions of phi which vary with position must be used.
% 2) The direction of micro-channels changes within the transition,
% entrance, and exit regions.  Roughly speaking the channels turn by 29
% degrees, and this serves as a simple estimation.  In reality the channels
% follow a more complex path, see InternalHeader_MicrochannelTransforms.pdf
% Transform functions, Gamma, will be used to rotate the permeability
% matrices within these regions.

% define core solid materials
% the core solid is stream 0 and we specify it first on face 3
sld1 = HHXTProperties(model,'Face',3,'Stream',0,...
                       'PhaseFraction',(1-phi_core),...
                 'ThermalConductivity',15.5,... %conductivity of 800H at 250 C in W/m-K
                         'MassDensity',8030,... %densitiy in kg/m3
                        'SpecificHeat',488.4);  %heat capacity in J/kg-K
                    
% the definition applied to Face 3 can be used as a basis for all other faces 
% by passing the HHXTMaterialAssignment object sld1.
% in the jump region
HHXTProperties(model,'Face',1,'Material',sld1,'PhaseFraction',@phi_S_jump);
HHXTProperties(model,'Face',10,'Material',sld1,'PhaseFraction',@phi_S_jump);
% in the cold entrance and exit regions
HHXTProperties(model,'Face',2,'Material',sld1,'PhaseFraction',@phi_SC_trans);
HHXTProperties(model,'Face',8,'Material',sld1,'PhaseFraction',@phi_SC_trans);
% in the hot engrance and exit regions
HHXTProperties(model,'Face',5,'Material',sld1,'PhaseFraction',@phi_SH_trans);
HHXTProperties(model,'Face',9,'Material',sld1,'PhaseFraction',@phi_SH_trans);
% in the transition region
HHXTProperties(model,'Face',4,'Material',sld1,'PhaseFraction',@phi_SHC_trans);
HHXTProperties(model,'Face',6,'Material',sld1,'PhaseFraction',@phi_SHC_trans);
% in the solid wall region the HX is only solid body material
HHXTProperties(model,'Face',7,'Material',sld1,'PhaseFraction',1);

% define cold stream (CO2) material
% the cold stream is stream 1 and we specify it first on face 3
% using real gas CO2 properties.
fldC = HHXTProperties(model,'Face',3,'Stream',1,...
                       'PhaseFraction',phi_C,...
                 'ThermalConductivity',@CO2_k_real,...       %thermal conductivity in W/m-K
                         'MassDensity',@CO2_rho_real,...     %densitiy in kg/m3
                        'SpecificHeat',@CO2_cp_real,...      %heat capacity in J/kg-K
                   'HydraulicDiameter',D_h,...      %hydraulic diameter in m
                           'Viscosity',@CO2_mu_real,...      %viscosity in kg/m-s
                             'Nusselt',Nu,...
                              'fDarcy',fmatr);
                          
% define hot stream (CO2) material
% the cold stream is stream 2 and we specify it first on face 3
% using real gas CO2 properties.
fldH = HHXTProperties(model,'Face',3,'Stream',2,...
                       'PhaseFraction',phi_H,...
                 'ThermalConductivity',@CO2_k_real,...       %thermal conductivity in W/m-K
                         'MassDensity',@CO2_rho_real,...     %densitiy in kg/m3
                        'SpecificHeat',@CO2_cp_real,...      %heat capacity in J/kg-K
                   'HydraulicDiameter',D_h,...      %hydraulic diameter in m
                           'Viscosity',@CO2_mu_real,...      %viscosity in kg/m-s
                             'Nusselt',Nu,...
                              'fDarcy',fmatr);
                          
% specifiy the jump regions
% the definitions applied to Face 3 can be copied to Faces 1 and 10 by
% passing the HHXTMaterialAssignment objects fldC and fldH 
HHXTProperties(model,'Face',1,'Material',fldC,'PhaseFraction',@phi_C_jump);
HHXTProperties(model,'Face',10,'Material',fldC,'PhaseFraction',@phi_C_jump);
HHXTProperties(model,'Face',1,'Material',fldH,'PhaseFraction',@phi_H_jump);
HHXTProperties(model,'Face',10,'Material',fldH,'PhaseFraction',@phi_H_jump);
                          
% specifiy the other regions based on GammaMode
% the definitions applied to Face 3 can be copied by
% passing the HHXTMaterialAssignment objects fldC and fldH                          
switch GammaMode

    case 'simple' % here we just apply a 29 degree rotation
    % determine the rotational transform matrices
    Gamma_pos29deg = GammafromTheta(29*(pi/180)); %for positive 29 deg
    Gamma_pos29deg = Gamma_pos29deg{1}; %have to extract from cell form
    Gamma_neg29deg = GammafromTheta(-29*(pi/180)); % for negative 29 deg
    Gamma_neg29deg = Gamma_neg29deg{1}; %have to extract from cell form
    % specify the transition regions on both the cold and hot sides
    HHXTProperties(model,'Face',4,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',Gamma_pos29deg);
    HHXTProperties(model,'Face',6,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',Gamma_pos29deg);
    HHXTProperties(model,'Face',4,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',Gamma_neg29deg);
    HHXTProperties(model,'Face',6,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',Gamma_neg29deg);
    % specify the entrance and exit regions for the cold side
    HHXTProperties(model,'Face',2,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',Gamma_pos29deg);
    HHXTProperties(model,'Face',8,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',Gamma_pos29deg);
    % specify the entrance and exit regions for the hot side
    HHXTProperties(model,'Face',5,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',Gamma_neg29deg);
    HHXTProperties(model,'Face',9,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',Gamma_neg29deg);

    case 'complex' % a more complex definition that depends on position
    % specify the transition regions on both the cold and hot sides
    HHXTProperties(model,'Face',4,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',@Gamma_C_4);
    HHXTProperties(model,'Face',6,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',@Gamma_C_6);
    HHXTProperties(model,'Face',4,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',@Gamma_H_4);
    HHXTProperties(model,'Face',6,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',@Gamma_H_6);
    % specify the entrance and exit regions for the cold side
    HHXTProperties(model,'Face',2,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',@Gamma_C_2);
    HHXTProperties(model,'Face',8,'Material',fldC,...
        'PhaseFraction',@phi_C_trans,'PermeabilityTransform',@Gamma_C_8);
    % specify the entrance and exit regions for the hot side
    HHXTProperties(model,'Face',5,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',@Gamma_H_5);
    HHXTProperties(model,'Face',9,'Material',fldH,...
        'PhaseFraction',@phi_H_trans,'PermeabilityTransform',@Gamma_H_9);

end

%% Apply Initial Conditions
% Apply initial conditions using the |setInitialConditions()| function.

% set initial conditions
% these are just guesses so the solver can start from somewhere
% we will guess average values
u0 = [0.25*(T_C_in+T_C_out+T_H_in+T_H_out),...
      0.5*(P_C_in+P_C_out),0.5*(T_C_in+T_C_out),...
      0.5*(P_H_in+P_H_out),0.5*(T_H_in+T_H_out)]';
setInitialConditions(model,u0);

%% Define Boundary Conditions
% Define boundary conditions using the |applyStreamBoundaryCondition()| function.  

% delete(model.BoundaryConditions);

% note that the inlets and outlets are each composed of two edges
C_edges_in = [36,37];
H_edges_out = [51,52];
H_edges_in = [43,44];
C_edges_out = [45,49];

switch flowBC
    case 'DP'
        %cold stream
        applyStreamBoundaryCondition(model,'Edge',C_edges_in,'Stream',1,'Direction','inlet',...
            'Pressure',P_C_in,'Temperature',T_C_in);
        applyStreamBoundaryCondition(model,'Edge',C_edges_out,'Stream',1,'Direction','outlet',...
            'Pressure',P_C_out);
        %hot stream
        applyStreamBoundaryCondition(model,'Edge',H_edges_in,'Stream',2,'Direction','inlet',...
            'Pressure',P_H_in,'Temperature',T_H_in);
        applyStreamBoundaryCondition(model,'Edge',H_edges_out,'Stream',2,'Direction','outlet',...
            'Pressure',P_H_out);
    case 'mdot'
        %cold stream
        applyStreamBoundaryCondition(model,'Edge',C_edges_in,'Stream',1,'Direction','inlet',...
            'Pressure',P_C_in,'Temperature',T_C_in,...
            'MassFlow',m_dot_C);
        applyStreamBoundaryCondition(model,'Edge',C_edges_out,'Stream',1,'Direction','outlet',...
            'MassFlow',m_dot_C);
        %hot stream
        applyStreamBoundaryCondition(model,'Edge',H_edges_in,'Stream',2,'Direction','inlet',...
            'Pressure',P_H_in,'Temperature',T_H_in,...
            'MassFlow',m_dot_H);
        applyStreamBoundaryCondition(model,'Edge',H_edges_out,'Stream',2,'Direction','outlet',...
            'MassFlow',m_dot_H);
end



%% Solve the steady state problem
% Solve the HHXT model using the |solveHHXTpde()| function.

if verbose == true
disp('starting steady state solution :...')    
end

%run the steady state solution, and time and display progress
if verbose == true
tic;
disp('    ...solving step')
end
results = solveHHXTpde(model); % solve the model
if verbose == true
solve_time = toc;
finish_date = datestr(now);
disp(['finished steady state solution at ',finish_date,' elapsed time ',num2str(solve_time),' sec'])
end

% save the model and results to file
save(filename,'model','results'); 
if verbose == true
disp(['    state saved to ',filename,'.mat']);
end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint
end

%% Display Results
% View the part of the results object.  The |results.Tables| provide an
% overview of the problem.

% display results
% results are given in results.Tables containing the following tables
% results.Tables{1}, a summary of the fluid stream conditions
% results.Tables{2}, solid body boundaries (none were specified in this
% model, so this is empty
% results.Tables{3}, stream 1 boundaries (cold stream, inlet and outlets)
% results.Tables{4}, stream 2 boundaries (hot stream, inlet and outlets)
if verbose == true
disp('Results summary:')
disp(results.Tables{1})

disp('')
disp('Cold Stream:')
disp(results.Tables{3})

disp('')
disp('Hot Stream:')
disp(results.Tables{4})

disp('')
disp('Number of Elements:')
disp(length(results.Mesh.Elements))
end

%% Plot Model Results
% Plot the results of the model using the |hhxt.PlotHHXT()} class.

if plotting == true
dir_plotting = [cd,'\figures'];
PlotResults(filename,dir_plotting,breakpoints);
end

%% Change in Volume Fraction Functions
% Functions defining the variation of the fluid phase fraction with the model. 
% Functional definitions must accept the inputs (pos,time)
% Where pos is a structure containing arrays of coordinates pos [x;y]
% And where time is an array of time in sec.
%
% 

function phi = phi_C_jump(pos,time)
%phi_C_jump - returns the volume fraction of the cold stream within the
%jump regions (faces 1 and 10)
    global phi_C
    phi = phi_C;

%     % the volume fraction of the cold channels decreases is halved across
%     % the jump region
%     x_7 = 0.304551080000000;    
%     x_8 = 0.308592220000000;
%     x = abs(pos(1,:));
%     val = interp1([x_8,x_7],[phi,0.5*phi],x);
%     val(x<x_8) = phi; val(x>x_7) = x_7;
%     phi = val; % cold channel volume fraction 
end

function phi = phi_H_jump(pos,time)
%phi_H_jump - returns the volume fraction of the cold stream within the
%jump regions (faces 1 and 10)
    global phi_H
    phi = phi_H;

%     % the volume fraction of the hot channels decreases is halved across
%     % the jump region
%     x_7 = 0.304551080000000;    
%     x_8 = 0.308592220000000;
%     x = abs(pos(1,:));
%     val = interp1([x_8,x_7],[phi,0.5*phi],x);
%     val(x<x_8) = phi; val(x>x_7) = x_7;
%     phi = val; % hot channel volume fraction 
end

function phi = phi_S_jump(pos,time)
%phi_S_jump - returns the volume fraction of the solid within the
%jump regions (faces 1 and 10)
    % solid channel volume fraction can be found from the 
    % hot and cold channel volume fractions
    phi = 1 - phi_H_jump(pos,time) - phi_C_jump(pos,time); 
end

function phi = phi_C_trans(pos,time)
%phi_C_trans -  returns the volume fraction of the cold stream within the
%transition, entrance, and exit regions (faces 2, 4, 6, and 8)
    global phi_C
    phi = 0.5*phi_C; % the volume fraction of cold channels is half that of the core
    
%     % the volume fraction of the hot channels incrreases by x2.2135 as we 
%     % get closer to the inlet/outlet
%     x_9 = 0.311134760000000;
%     x_10 = 0.342752680000000;
%     x = abs(pos(1,:));
%     val = interp1([x_9,x_10],[phi,2.2135*phi],x);
%     val(x<=x_9) = phi; val(x>=x_10) = 2.2135*phi;
%     phi = val; % cold channel volume fraction 
end

function phi = phi_H_trans(pos,time)
%phi_H_trans -  returns the volume fraction of the hot stream within the
%entrance and exit regions (faces 4, 5, 6, and 9)
    global phi_H
    phi = 0.5*phi_H; % the volume fraction of hot channels is half that of the core
    
%     % the volume fraction of the hot channels incrreases by x2.2135 as we 
%     % get closer to the inlet/outlet
%     x_9 = 0.311134760000000;
%     x_10 = 0.342752680000000;
%     x = abs(pos(1,:));
%     val = interp1([x_9,x_10],[phi,2.2135*phi],x);
%     val(x<=x_9) = phi; val(x>=x_10) = 2.2135*phi;
%     phi = val; % hot channel volume fraction 
end

function phi = phi_SHC_trans(pos,time)
%phi_SHC_trans -  returns the volume fraction of the solid within the
%transition regions (faces 4 and 6)

    % solid channel volume fraction can be found from the 
    % hot and cold channel volume fractions
    phi = 1 - phi_H_trans(pos,time) - phi_C_trans(pos,time); 
end

function phi = phi_SC_trans(pos,time)
%phi_SC_trans -  returns the volume fraction of the solid within the cold
%entrance and exit regions (faces 2 and 8)
    
    % solid channel volume fraction can be found from the 
    % cold channel volume fraction 
    phi = 1 - phi_C_trans(pos,time); 
end

function phi = phi_SH_trans(pos,time)
%phi_SH_trans -  returns the volume fraction of the solid within the hot
%entrance and exit regions (faces 5 and 9)
    
    % solid channel volume fraction can be found from the 
    % hot channel volume fraction 
    phi = 1 - phi_H_trans(pos,time); 
end

%% Region Transform Functions
% Transformation of the micro-channel direction can depend on the model and
% positions.
% Functional definitions must accept the inputs (model,pos)
% Where model is the HHXT model and pos is a structure containing arrays 
% of coordinates pos [x;y].  The model is passed so that mesh data or
% geometric features can be used (model.Geometry and model.Mesh).
%
% Which transform functions are called depends on the choice of GammaMode
% if GammaMode = 'complex', a multiple rotation form is applied that
%       accurately follows the actual layout of the etched channels.  This
%       utilizes the Gamma_C_* and Gamma_H_* functions.
% see InternalHeader_MicrochannelTransforms.pdf for a diagram of internal
% header micro-channels and the related transformation function variables.

function Gamma = Gamma_C_4(model,pos)
%Gamma_C_4 - returns the rotational transform matrices for the cold side
%within face 4.
    theta = rotate1(pos(1,:),pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_H_4(model,pos)
%Gamma_H_4 - returns the rotational transform matrices for the hot side
%within face 4.
    theta = -rotate1(pos(1,:),-pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_C_6(model,pos)
%Gamma_C_6 - returns the rotational transform matrices for the cold side
%within face 6.
    theta = rotate1(-pos(1,:),-pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_H_6(model,pos)
%Gamma_H_6 - returns the rotational transform matrices for the hot side
%within face 6.
    theta = -rotate1(-pos(1,:),pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_C_2(model,pos)
%Gamma_C_2 - returns the rotational transform matrices for the cold side
%within face 2.
    theta = rotate2(pos(1,:),pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_H_9(model,pos)
%Gamma_H_9 - returns the rotational transform matrices for the hot side
%within face 9.
    theta = -rotate2(pos(1,:),-pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_C_8(model,pos)
%Gamma_C_8 - returns the rotational transform matrices for the cold side
%within face 8.
    theta = rotate2(-pos(1,:),-pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = Gamma_H_5(model,pos)
%Gamma_H_5 - returns the rotational transform matrices for the hot side
%within face 5.
    theta = -rotate2(-pos(1,:),pos(2,:)); % find angle of rotation based on x y position
    Gamma = GammafromTheta(theta); % get the transform matrices for theta rotation about z
end

function Gamma = GammafromTheta(theta)
%GammafromTheta - returns a cell array containing rotational transform
%matrices Gamma which have been rotate by theta radians about the z axis
    c = cos(theta); s = sin(theta);
    Gamma = cell(1,length(theta));
    for k = 1:length(theta)
    Gamma{k} = [c(k),-s(k);s(k),c(k)];
    end
end

function theta = rotate1(x,y)
%rotate1 - returns the angle of rotation to be applied to positions of x, y
% with respect to the rotational center c_2 = [-.3665,-0.0429]

	normal = [1,0]; % the normal of the unrotated space, in x direction, i.e. global norm
    
    % define the center of rotation
    % this is the intersection of edge 6 and edge 10
    x_2 = -0.366450880000000;
    y_2 = -0.042892980000000;
    centers = ones(length(x),1)*[x_2,y_2]; 
    
    % determine the vector norm of the vector pointing from the centers to
    % the x,y position
    r_vec = [x',y'] - centers; %the vector pointing from the centers to the x, y posistions
    r_mag = ones(length(x),1); %the magnitude of r_vec
    for k = 1:length(x)
        r_mag(k,1) = norm(r_vec(k,:));
    end
    
    % use definition of dot product of r_vec and normal to find theta
    theta = acos(r_vec*normal'./r_mag); 
    
    % set theta to 0 for any positions to the right of x_6
    x_6 = -0.311134760000000;
    theta(x>x_6) = 0;    
end

function theta = rotate2(x,y)
%rotate2 - returns the angle of rotation to be applied to positions of x, y
% with respect to either rotational center c_3 = [-.3556,-0.0176] or
% c_1 = [-.3799,-.0429] depending on 

    normal = [1,0]; % the normal of the unrotated space, in x direction, i.e. global norm
     
    % define the center of rotation
    % this is the intersection of edge 17 and edge 18
    x_3 = -0.355640640000000;
    y_3 = -0.0175971200000000;
    centers = ones(length(x),1)*[x_3,y_3];
    
    % define the center of rotation
    % this is the intersection of edge 10 and edge 18
    x_1 = -0.379859540000000;
    y_1 = -0.042892980000000;
    
    % define coordinates of a line between vertex 4 and vertex 9
    x_4 = -0.342752680000000;
    y_4 = y_1;
    x_5 = -0.338792820000000;
    y_5 = 0;
    % from the y positions, determine which x positions are on the line
    x_line = ( (y-y_1)/(-y_1) )*(x_5-x_4)+x_4; 
    
    % if the x y position is to the right of the line defined by x_line
    % then the center of rotation is about x_1 y_1
    centers(x>=x_line,:) = ones(sum(x>=x_line),1)*[x_1,y_1];
    
    % determine the vector norm of the vector pointing from the centers to
    % the x,y position
    r_vec = [x',y'] - centers; %the vector pointing from the centers to the x, y posistions
    r_mag = ones(length(x),1); %the magnitude of r_vec
    for k = 1:length(x)
        r_mag(k,1) = norm(r_vec(k,:));
    end
    
    % use definition of dot product of r_vec and normal to find theta
    theta = acos(r_vec*normal'./r_mag); 
    
%     theta(y'<centers(:,2)) = -theta(y'<centers(:,2));

end

%% RealGas Fluid Properties: Varying Pressure and Temperature
% Properties of the fluids, streams 1 & 2, can be dependent on position,
% time, temperature, and pressure.  
% Functional definitions must accept the inputs (pos,time,temp,press)
% Where pos is a structure containing arrays of coordinates pos [x;y].
% And where time, temp, and press are arrays of time in sec,
% temperature in C, and pressure in Pa.

function rho = CO2_rho_real(pos,time,temp,press)
    %for 7.5 < P < 32.5 MPa, 20 < T < 700 C
    T = [20;33.6000000000000;47.2000000000000;60.8000000000000;74.4000000000000;88;101.600000000000;115.200000000000;128.800000000000;142.400000000000;156;169.600000000000;183.200000000000;196.800000000000;210.400000000000;224;237.600000000000;251.200000000000;264.800000000000;278.400000000000;292;305.600000000000;319.200000000000;332.800000000000;346.400000000000;360;373.600000000000;387.200000000000;400.800000000000;414.400000000000;428;441.600000000000;455.200000000000;468.800000000000;482.400000000000;496;509.600000000000;523.200000000000;536.800000000000;550.400000000000;564;577.600000000000;591.200000000000;604.800000000000;618.400000000000;632;645.600000000000;659.200000000000;672.800000000000;686.400000000000;700];
    P = [7.50000000000000,8.50000000000000,9,9.50000000000000,10,10.5000000000000,11,12,13,13.8000000000000,14,15,17,19,21,23,25,27.5000000000000,30,32.5000000000000];
    P = P*1e6;
    [T,P] = meshgrid(T,P);
    y = [818.700000000000,835.800000000000,843.200000000000,850,856.300000000000,862.200000000000,867.800000000000,878.100000000000,887.400000000000,894.400000000000,896,904,918.300000000000,931.100000000000,942.700000000000,953.300000000000,963,974.300000000000,984.700000000000,994.400000000000;296.200000000000,655.600000000000,689.900000000000,713.100000000000,731,745.800000000000,758.400000000000,779.400000000000,796.700000000000,808.600000000000,811.400000000000,824.300000000000,846.200000000000,864.600000000000,880.500000000000,894.500000000000,907.200000000000,921.500000000000,934.400000000000,946.200000000000;201.900000000000,265.200000000000,309.700000000000,368.100000000000,439.800000000000,508.100000000000,560.300000000000,627.200000000000,669.400000000000,694.500000000000,699.900000000000,723.900000000000,760.600000000000,788.500000000000,811.200000000000,830.400000000000,847.100000000000,865.400000000000,881.600000000000,896.100000000000;171.400000000000,210.200000000000,232.700000000000,257.700000000000,285.500000000000,316.600000000000,350.800000000000,425,495.400000000000,542,552.200000000000,596,658.200000000000,701.300000000000,734,760.500000000000,782.700000000000,806.200000000000,826.400000000000,844.200000000000;152.800000000000,182.700000000000,199.100000000000,216.500000000000,235,254.700000000000,275.700000000000,321.200000000000,370.400000000000,410.600000000000,420.500000000000,468.100000000000,547.600000000000,606.500000000000,650.900000000000,686,714.700000000000,744.400000000000,769.400000000000,790.900000000000;139.500000000000,164.600000000000,177.900000000000,191.800000000000,206.300000000000,221.400000000000,237.100000000000,270.400000000000,306,335.800000000000,343.300000000000,381.300000000000,454.300000000000,517.500000000000,569.300000000000,611.300000000000,646,681.900000000000,711.700000000000,737.100000000000;129.300000000000,151.100000000000,162.600000000000,174.400000000000,186.600000000000,199.100000000000,212,239,267.300000000000,290.800000000000,296.800000000000,327.100000000000,388.200000000000,446.100000000000,498,542.900000000000,581.200000000000,621.500000000000,655.400000000000,684.300000000000;121,140.600000000000,150.700000000000,161.100000000000,171.700000000000,182.600000000000,193.800000000000,216.800000000000,240.600000000000,260.300000000000,265.300000000000,290.600000000000,342.100000000000,392.900000000000,440.900000000000,484.700000000000,523.700000000000,566.200000000000,602.700000000000,634.200000000000;114,131.900000000000,141.100000000000,150.400000000000,160,169.700000000000,179.600000000000,199.900000000000,220.700000000000,237.800000000000,242.200000000000,264,308.600000000000,353.200000000000,396.500000000000,437.400000000000,475.100000000000,517.600000000000,555,587.900000000000;108.100000000000,124.600000000000,133,141.600000000000,150.300000000000,159.100000000000,168.100000000000,186.300000000000,205,220.300000000000,224.100000000000,243.500000000000,283,322.700000000000,361.700000000000,399.300000000000,434.900000000000,475.900000000000,513,546.300000000000;102.900000000000,118.300000000000,126.100000000000,134.100000000000,142.100000000000,150.200000000000,158.500000000000,175.200000000000,192.200000000000,206,209.500000000000,227.100000000000,262.700000000000,298.500000000000,333.900000000000,368.500000000000,401.600000000000,440.600000000000,476.600000000000,509.500000000000;98.3400000000000,112.800000000000,120.100000000000,127.500000000000,135,142.600000000000,150.200000000000,165.700000000000,181.400000000000,194.100000000000,197.300000000000,213.400000000000,246,278.700000000000,311.200000000000,343.100000000000,373.900000000000,410.700000000000,445.200000000000,477.200000000000;94.2500000000000,107.900000000000,114.800000000000,121.800000000000,128.800000000000,135.900000000000,143.100000000000,157.500000000000,172.100000000000,183.900000000000,186.900000000000,201.800000000000,231.900000000000,262.200000000000,292.200000000000,321.800000000000,350.600000000000,385.300000000000,418.100000000000,449;90.5700000000000,103.500000000000,110.100000000000,116.700000000000,123.300000000000,130,136.800000000000,150.300000000000,164.100000000000,175.100000000000,177.900000000000,191.800000000000,219.900000000000,248,276,303.700000000000,330.700000000000,363.400000000000,394.600000000000,424.200000000000;87.2300000000000,99.5900000000000,105.800000000000,112.100000000000,118.400000000000,124.800000000000,131.100000000000,144,156.900000000000,167.300000000000,170,183.100000000000,209.400000000000,235.800000000000,262.100000000000,288,313.500000000000,344.400000000000,374.100000000000,402.400000000000;84.1800000000000,96,102,107.900000000000,114,120,126.100000000000,138.300000000000,150.500000000000,160.400000000000,162.900000000000,175.300000000000,200.200000000000,225.100000000000,249.800000000000,274.300000000000,298.400000000000,327.700000000000,356,383.200000000000;81.3700000000000,92.7100000000000,98.4100000000000,104.100000000000,109.900000000000,115.700000000000,121.500000000000,133.100000000000,144.800000000000,154.200000000000,156.500000000000,168.300000000000,191.900000000000,215.500000000000,239,262.200000000000,285.100000000000,313,340,366;78.7700000000000,89.6800000000000,95.1600000000000,100.700000000000,106.200000000000,111.700000000000,117.300000000000,128.400000000000,139.600000000000,148.500000000000,150.800000000000,162,184.500000000000,207,229.300000000000,251.400000000000,273.200000000000,299.800000000000,325.700000000000,350.700000000000;76.3700000000000,86.8800000000000,92.1500000000000,97.4500000000000,102.800000000000,108.100000000000,113.400000000000,124.100000000000,134.800000000000,143.400000000000,145.500000000000,156.300000000000,177.800000000000,199.300000000000,220.600000000000,241.700000000000,262.500000000000,288,312.800000000000,336.800000000000;74.1200000000000,84.2700000000000,89.3700000000000,94.4700000000000,99.5800000000000,104.700000000000,109.800000000000,120.100000000000,130.400000000000,138.700000000000,140.700000000000,151.100000000000,171.700000000000,192.300000000000,212.700000000000,232.900000000000,252.800000000000,277.300000000000,301.100000000000,324.200000000000;72.0300000000000,81.8500000000000,86.7700000000000,91.7000000000000,96.6400000000000,101.600000000000,106.500000000000,116.400000000000,126.400000000000,134.300000000000,136.300000000000,146.200000000000,166.100000000000,185.800000000000,205.400000000000,224.800000000000,244,267.500000000000,290.500000000000,312.800000000000;70.0600000000000,79.5800000000000,84.3400000000000,89.1100000000000,93.8900000000000,98.6700000000000,103.500000000000,113,122.600000000000,130.300000000000,132.200000000000,141.800000000000,160.900000000000,179.900000000000,198.800000000000,217.500000000000,235.900000000000,258.600000000000,280.700000000000,302.300000000000;68.2200000000000,77.4500000000000,82.0700000000000,86.6900000000000,91.3200000000000,95.9500000000000,100.600000000000,109.800000000000,119.100000000000,126.500000000000,128.400000000000,137.600000000000,156.100000000000,174.500000000000,192.700000000000,210.700000000000,228.500000000000,250.300000000000,271.700000000000,292.600000000000;66.4800000000000,75.4400000000000,79.9300000000000,84.4200000000000,88.9100000000000,93.4000000000000,97.8900000000000,106.900000000000,115.900000000000,123,124.800000000000,133.800000000000,151.600000000000,169.400000000000,187,204.400000000000,221.600000000000,242.700000000000,263.400000000000,283.700000000000;64.8400000000000,73.5500000000000,77.9100000000000,82.2700000000000,86.6400000000000,91,95.3600000000000,104.100000000000,112.800000000000,119.700000000000,121.500000000000,130.200000000000,147.500000000000,164.700000000000,181.700000000000,198.600000000000,215.200000000000,235.700000000000,255.800000000000,275.400000000000;63.2800000000000,71.7700000000000,76.0100000000000,80.2500000000000,84.4900000000000,88.7300000000000,92.9700000000000,101.400000000000,109.900000000000,116.700000000000,118.400000000000,126.800000000000,143.600000000000,160.200000000000,176.800000000000,193.100000000000,209.200000000000,229.100000000000,248.600000000000,267.600000000000;61.8100000000000,70.0700000000000,74.2100000000000,78.3400000000000,82.4700000000000,86.6000000000000,90.7200000000000,98.9700000000000,107.200000000000,113.800000000000,115.400000000000,123.600000000000,139.900000000000,156.100000000000,172.100000000000,188,203.700000000000,223,241.900000000000,260.400000000000;60.4100000000000,68.4700000000000,72.5000000000000,76.5300000000000,80.5500000000000,84.5700000000000,88.5900000000000,96.6200000000000,104.600000000000,111,112.600000000000,120.600000000000,136.500000000000,152.200000000000,167.800000000000,183.200000000000,198.500000000000,217.200000000000,235.600000000000,253.700000000000;59.0700000000000,66.9400000000000,70.8800000000000,74.8000000000000,78.7300000000000,82.6500000000000,86.5700000000000,94.4000000000000,102.200000000000,108.400000000000,110,117.800000000000,133.200000000000,148.500000000000,163.700000000000,178.700000000000,193.500000000000,211.800000000000,229.700000000000,247.300000000000;57.8000000000000,65.4900000000000,69.3300000000000,73.1700000000000,77,80.8300000000000,84.6500000000000,92.2900000000000,99.9000000000000,106,107.500000000000,115.100000000000,130.100000000000,145,159.800000000000,174.500000000000,188.900000000000,206.700000000000,224.200000000000,241.300000000000;56.5900000000000,64.1100000000000,67.8600000000000,71.6100000000000,75.3500000000000,79.0900000000000,82.8300000000000,90.2800000000000,97.7200000000000,103.600000000000,105.100000000000,112.500000000000,127.200000000000,141.800000000000,156.200000000000,170.400000000000,184.500000000000,201.900000000000,219,235.700000000000;55.4300000000000,62.7800000000000,66.4500000000000,70.1200000000000,73.7800000000000,77.4300000000000,81.0800000000000,88.3700000000000,95.6300000000000,101.400000000000,102.900000000000,110.100000000000,124.400000000000,138.600000000000,152.700000000000,166.600000000000,180.400000000000,197.400000000000,214,230.300000000000;54.3300000000000,61.5200000000000,65.1100000000000,68.7000000000000,72.2800000000000,75.8500000000000,79.4200000000000,86.5500000000000,93.6500000000000,99.3100000000000,100.700000000000,107.800000000000,121.800000000000,135.700000000000,149.400000000000,163,176.500000000000,193,209.300000000000,225.300000000000;53.2700000000000,60.3100000000000,63.8300000000000,67.3400000000000,70.8400000000000,74.3400000000000,77.8400000000000,84.8100000000000,91.7600000000000,97.3000000000000,98.6800000000000,105.600000000000,119.300000000000,132.900000000000,146.300000000000,159.600000000000,172.700000000000,188.900000000000,204.800000000000,220.500000000000;52.2500000000000,59.1500000000000,62.5900000000000,66.0300000000000,69.4700000000000,72.8900000000000,76.3200000000000,83.1400000000000,89.9400000000000,95.3700000000000,96.7200000000000,103.500000000000,116.900000000000,130.200000000000,143.300000000000,156.300000000000,169.200000000000,185,200.600000000000,215.900000000000;51.2700000000000,58.0400000000000,61.4100000000000,64.7800000000000,68.1500000000000,71.5100000000000,74.8600000000000,81.5500000000000,88.2100000000000,93.5200000000000,94.8500000000000,101.500000000000,114.600000000000,127.600000000000,140.500000000000,153.200000000000,165.800000000000,181.300000000000,196.600000000000,211.600000000000;50.3300000000000,56.9700000000000,60.2800000000000,63.5800000000000,66.8800000000000,70.1800000000000,73.4600000000000,80.0200000000000,86.5500000000000,91.7600000000000,93.0500000000000,99.5300000000000,112.400000000000,125.100000000000,137.800000000000,150.200000000000,162.600000000000,177.800000000000,192.700000000000,207.400000000000;49.4300000000000,55.9400000000000,59.1900000000000,62.4300000000000,65.6700000000000,68.9000000000000,72.1200000000000,78.5500000000000,84.9600000000000,90.0600000000000,91.3300000000000,97.6800000000000,110.300000000000,122.800000000000,135.100000000000,147.400000000000,159.500000000000,174.400000000000,189,203.400000000000;48.5600000000000,54.9500000000000,58.1400000000000,61.3200000000000,64.5000000000000,67.6700000000000,70.8300000000000,77.1400000000000,83.4300000000000,88.4300000000000,89.6800000000000,95.9100000000000,108.300000000000,120.500000000000,132.700000000000,144.600000000000,156.500000000000,171.100000000000,185.500000000000,199.600000000000;47.7200000000000,54,57.1300000000000,60.2600000000000,63.3800000000000,66.4900000000000,69.5900000000000,75.7900000000000,81.9500000000000,86.8700000000000,88.1000000000000,94.2100000000000,106.300000000000,118.400000000000,130.300000000000,142,153.700000000000,168,182.100000000000,196;46.9200000000000,53.0900000000000,56.1600000000000,59.2300000000000,62.2900000000000,65.3500000000000,68.4000000000000,74.4800000000000,80.5400000000000,85.3600000000000,86.5700000000000,92.5700000000000,104.500000000000,116.300000000000,128,139.500000000000,150.900000000000,165,178.900000000000,192.500000000000;46.1400000000000,52.2000000000000,55.2200000000000,58.2400000000000,61.2500000000000,64.2500000000000,67.2500000000000,73.2300000000000,79.1700000000000,83.9200000000000,85.1000000000000,90.9900000000000,102.700000000000,114.300000000000,125.800000000000,137.100000000000,148.300000000000,162.100000000000,175.800000000000,189.200000000000;45.3900000000000,51.3500000000000,54.3200000000000,57.2800000000000,60.2400000000000,63.1900000000000,66.1400000000000,72.0100000000000,77.8600000000000,82.5200000000000,83.6800000000000,89.4700000000000,101,112.400000000000,123.600000000000,134.800000000000,145.800000000000,159.400000000000,172.800000000000,185.900000000000;44.6600000000000,50.5200000000000,53.4400000000000,56.3600000000000,59.2700000000000,62.1700000000000,65.0700000000000,70.8400000000000,76.5900000000000,81.1700000000000,82.3100000000000,88.0100000000000,99.3200000000000,110.500000000000,121.600000000000,132.500000000000,143.400000000000,156.700000000000,169.900000000000,182.800000000000;43.9600000000000,49.7200000000000,52.6000000000000,55.4700000000000,58.3300000000000,61.1800000000000,64.0300000000000,69.7100000000000,75.3700000000000,79.8700000000000,80.9900000000000,86.6000000000000,97.7100000000000,108.700000000000,119.600000000000,130.400000000000,141,154.200000000000,167.100000000000,179.900000000000;43.2800000000000,48.9500000000000,51.7800000000000,54.6000000000000,57.4200000000000,60.2300000000000,63.0300000000000,68.6200000000000,74.1800000000000,78.6200000000000,79.7200000000000,85.2300000000000,96.1700000000000,107,117.700000000000,128.300000000000,138.800000000000,151.700000000000,164.400000000000,177;42.6200000000000,48.2100000000000,50.9900000000000,53.7700000000000,56.5400000000000,59.3100000000000,62.0700000000000,67.5700000000000,73.0400000000000,77.4000000000000,78.4900000000000,83.9100000000000,94.6700000000000,105.300000000000,115.900000000000,126.300000000000,136.600000000000,149.300000000000,161.900000000000,174.200000000000;41.9800000000000,47.4800000000000,50.2200000000000,52.9600000000000,55.6900000000000,58.4100000000000,61.1300000000000,66.5400000000000,71.9300000000000,76.2300000000000,77.3000000000000,82.6300000000000,93.2300000000000,103.700000000000,114.100000000000,124.400000000000,134.500000000000,147,159.400000000000,171.500000000000;41.3600000000000,46.7800000000000,49.4800000000000,52.1800000000000,54.8600000000000,57.5500000000000,60.2200000000000,65.5500000000000,70.8600000000000,75.0900000000000,76.1400000000000,81.4000000000000,91.8300000000000,102.200000000000,112.400000000000,122.500000000000,132.500000000000,144.800000000000,157,168.900000000000;40.7600000000000,46.1000000000000,48.7600000000000,51.4200000000000,54.0600000000000,56.7100000000000,59.3400000000000,64.6000000000000,69.8200000000000,73.9900000000000,75.0300000000000,80.2000000000000,90.4800000000000,100.700000000000,110.700000000000,120.700000000000,130.500000000000,142.700000000000,154.600000000000,166.400000000000;40.1800000000000,45.4400000000000,48.0700000000000,50.6800000000000,53.2900000000000,55.8900000000000,58.4900000000000,63.6700000000000,68.8200000000000,72.9200000000000,73.9400000000000,79.0400000000000,89.1700000000000,99.1900000000000,109.100000000000,118.900000000000,128.600000000000,140.600000000000,152.400000000000,164];
    y = y';
    val = interpinbounds(T,P,y,temp,press);
    rho = val;
end

function mu = CO2_mu_real(pos,time,temp,press)
    %for 7.5 < P < 32.5 MPa, 20 < T < 700 C
    T = [20;33.6000000000000;47.2000000000000;60.8000000000000;74.4000000000000;88;101.600000000000;115.200000000000;128.800000000000;142.400000000000;156;169.600000000000;183.200000000000;196.800000000000;210.400000000000;224;237.600000000000;251.200000000000;264.800000000000;278.400000000000;292;305.600000000000;319.200000000000;332.800000000000;346.400000000000;360;373.600000000000;387.200000000000;400.800000000000;414.400000000000;428;441.600000000000;455.200000000000;468.800000000000;482.400000000000;496;509.600000000000;523.200000000000;536.800000000000;550.400000000000;564;577.600000000000;591.200000000000;604.800000000000;618.400000000000;632;645.600000000000;659.200000000000;672.800000000000;686.400000000000;700];
    P = [7.50000000000000,8.50000000000000,9,9.50000000000000,10,10.5000000000000,11,12,13,13.8000000000000,14,15,17,19,21,23,25,27.5000000000000,30,32.5000000000000];
    P = P*1e6;
    [T,P] = meshgrid(T,P);
    y = [7.40200000000000e-05,7.72900000000000e-05,7.87700000000000e-05,8.01600000000000e-05,8.14900000000000e-05,8.27600000000000e-05,8.39700000000000e-05,8.62800000000000e-05,8.84600000000000e-05,9.01100000000000e-05,9.05200000000000e-05,9.24900000000000e-05,9.62000000000000e-05,9.96700000000000e-05,0.000103000000000000,0.000106100000000000,0.000109100000000000,0.000112800000000000,0.000116300000000000,0.000119600000000000;2.28700000000000e-05,5.05600000000000e-05,5.46700000000000e-05,5.76800000000000e-05,6.01400000000000e-05,6.22600000000000e-05,6.41400000000000e-05,6.74400000000000e-05,7.03100000000000e-05,7.23900000000000e-05,7.28900000000000e-05,7.52500000000000e-05,7.95000000000000e-05,8.33100000000000e-05,8.68100000000000e-05,9.00800000000000e-05,9.31800000000000e-05,9.68400000000000e-05,0.000100300000000000,0.000103700000000000;1.97500000000000e-05,2.21400000000000e-05,2.41700000000000e-05,2.72800000000000e-05,3.18100000000000e-05,3.69200000000000e-05,4.13900000000000e-05,4.79700000000000e-05,5.27000000000000e-05,5.57500000000000e-05,5.64500000000000e-05,5.96100000000000e-05,6.48800000000000e-05,6.92800000000000e-05,7.31500000000000e-05,7.66500000000000e-05,7.98900000000000e-05,8.36500000000000e-05,8.71700000000000e-05,9.05100000000000e-05;1.94500000000000e-05,2.06800000000000e-05,2.14900000000000e-05,2.24700000000000e-05,2.36800000000000e-05,2.51600000000000e-05,2.69500000000000e-05,3.14500000000000e-05,3.65300000000000e-05,4.03800000000000e-05,4.12800000000000e-05,4.53800000000000e-05,5.19500000000000e-05,5.71400000000000e-05,6.14900000000000e-05,6.53000000000000e-05,6.87400000000000e-05,7.26400000000000e-05,7.62300000000000e-05,7.95900000000000e-05;1.95800000000000e-05,2.04300000000000e-05,2.09400000000000e-05,2.15400000000000e-05,2.22100000000000e-05,2.29900000000000e-05,2.38800000000000e-05,2.60300000000000e-05,2.86900000000000e-05,3.11400000000000e-05,3.17900000000000e-05,3.51000000000000e-05,4.14900000000000e-05,4.70100000000000e-05,5.17000000000000e-05,5.57700000000000e-05,5.93900000000000e-05,6.34400000000000e-05,6.71100000000000e-05,7.05000000000000e-05;1.98800000000000e-05,2.05300000000000e-05,2.09200000000000e-05,2.13400000000000e-05,2.18200000000000e-05,2.23400000000000e-05,2.29300000000000e-05,2.42800000000000e-05,2.59100000000000e-05,2.74100000000000e-05,2.78100000000000e-05,2.99600000000000e-05,3.47200000000000e-05,3.95500000000000e-05,4.40400000000000e-05,4.80800000000000e-05,5.17300000000000e-05,5.58300000000000e-05,5.95400000000000e-05,6.29300000000000e-05;2.02700000000000e-05,2.08000000000000e-05,2.11100000000000e-05,2.14400000000000e-05,2.18100000000000e-05,2.22100000000000e-05,2.26400000000000e-05,2.36300000000000e-05,2.47800000000000e-05,2.58200000000000e-05,2.61000000000000e-05,2.75900000000000e-05,3.10000000000000e-05,3.47700000000000e-05,3.86100000000000e-05,4.22900000000000e-05,4.57400000000000e-05,4.97200000000000e-05,5.33500000000000e-05,5.66900000000000e-05;2.07100000000000e-05,2.11600000000000e-05,2.14100000000000e-05,2.16900000000000e-05,2.19900000000000e-05,2.23100000000000e-05,2.26600000000000e-05,2.34300000000000e-05,2.43100000000000e-05,2.51000000000000e-05,2.53100000000000e-05,2.64300000000000e-05,2.90000000000000e-05,3.19100000000000e-05,3.50300000000000e-05,3.81900000000000e-05,4.12800000000000e-05,4.49700000000000e-05,4.84100000000000e-05,5.16200000000000e-05;2.11800000000000e-05,2.15700000000000e-05,2.17900000000000e-05,2.20200000000000e-05,2.22700000000000e-05,2.25400000000000e-05,2.28300000000000e-05,2.34700000000000e-05,2.41800000000000e-05,2.48200000000000e-05,2.49900000000000e-05,2.58700000000000e-05,2.79000000000000e-05,3.02100000000000e-05,3.27400000000000e-05,3.54000000000000e-05,3.80800000000000e-05,4.13800000000000e-05,4.45600000000000e-05,4.75800000000000e-05;2.16700000000000e-05,2.20100000000000e-05,2.22000000000000e-05,2.24100000000000e-05,2.26200000000000e-05,2.28600000000000e-05,2.31000000000000e-05,2.36400000000000e-05,2.42500000000000e-05,2.47800000000000e-05,2.49200000000000e-05,2.56500000000000e-05,2.73100000000000e-05,2.92000000000000e-05,3.12900000000000e-05,3.35200000000000e-05,3.58300000000000e-05,3.87400000000000e-05,4.16100000000000e-05,4.44000000000000e-05;2.21700000000000e-05,2.24800000000000e-05,2.26500000000000e-05,2.28300000000000e-05,2.30200000000000e-05,2.32200000000000e-05,2.34400000000000e-05,2.39100000000000e-05,2.44300000000000e-05,2.48800000000000e-05,2.50000000000000e-05,2.56200000000000e-05,2.70200000000000e-05,2.86200000000000e-05,3.03800000000000e-05,3.22800000000000e-05,3.42700000000000e-05,3.68200000000000e-05,3.93900000000000e-05,4.19200000000000e-05;2.26800000000000e-05,2.29600000000000e-05,2.31100000000000e-05,2.32700000000000e-05,2.34400000000000e-05,2.36300000000000e-05,2.38200000000000e-05,2.42300000000000e-05,2.46900000000000e-05,2.50800000000000e-05,2.51900000000000e-05,2.57300000000000e-05,2.69400000000000e-05,2.83100000000000e-05,2.98300000000000e-05,3.14600000000000e-05,3.31900000000000e-05,3.54300000000000e-05,3.77300000000000e-05,4.00200000000000e-05;2.31900000000000e-05,2.34500000000000e-05,2.35900000000000e-05,2.37300000000000e-05,2.38900000000000e-05,2.40500000000000e-05,2.42300000000000e-05,2.46000000000000e-05,2.50000000000000e-05,2.53500000000000e-05,2.54500000000000e-05,2.59200000000000e-05,2.69900000000000e-05,2.81900000000000e-05,2.95100000000000e-05,3.09400000000000e-05,3.24600000000000e-05,3.44400000000000e-05,3.64900000000000e-05,3.85700000000000e-05;2.37100000000000e-05,2.39400000000000e-05,2.40700000000000e-05,2.42100000000000e-05,2.43500000000000e-05,2.45000000000000e-05,2.46500000000000e-05,2.49900000000000e-05,2.53600000000000e-05,2.56700000000000e-05,2.57500000000000e-05,2.61800000000000e-05,2.71300000000000e-05,2.81900000000000e-05,2.93600000000000e-05,3.06300000000000e-05,3.19800000000000e-05,3.37500000000000e-05,3.55900000000000e-05,3.74700000000000e-05;2.42300000000000e-05,2.44500000000000e-05,2.45600000000000e-05,2.46900000000000e-05,2.48200000000000e-05,2.49500000000000e-05,2.51000000000000e-05,2.54000000000000e-05,2.57400000000000e-05,2.60200000000000e-05,2.61000000000000e-05,2.64800000000000e-05,2.73400000000000e-05,2.82900000000000e-05,2.93400000000000e-05,3.04700000000000e-05,3.16800000000000e-05,3.32700000000000e-05,3.49400000000000e-05,3.66500000000000e-05;2.47500000000000e-05,2.49500000000000e-05,2.50600000000000e-05,2.51700000000000e-05,2.52900000000000e-05,2.54200000000000e-05,2.55500000000000e-05,2.58300000000000e-05,2.61400000000000e-05,2.64000000000000e-05,2.64700000000000e-05,2.68200000000000e-05,2.75900000000000e-05,2.84600000000000e-05,2.94100000000000e-05,3.04300000000000e-05,3.15200000000000e-05,3.29700000000000e-05,3.44800000000000e-05,3.60400000000000e-05;2.52600000000000e-05,2.54500000000000e-05,2.55500000000000e-05,2.56600000000000e-05,2.57700000000000e-05,2.58900000000000e-05,2.60100000000000e-05,2.62700000000000e-05,2.65600000000000e-05,2.68000000000000e-05,2.68600000000000e-05,2.71800000000000e-05,2.78900000000000e-05,2.86800000000000e-05,2.95400000000000e-05,3.04700000000000e-05,3.14700000000000e-05,3.27900000000000e-05,3.41700000000000e-05,3.56100000000000e-05;2.57800000000000e-05,2.59600000000000e-05,2.60500000000000e-05,2.61500000000000e-05,2.62600000000000e-05,2.63700000000000e-05,2.64800000000000e-05,2.67200000000000e-05,2.69800000000000e-05,2.72100000000000e-05,2.72600000000000e-05,2.75600000000000e-05,2.82100000000000e-05,2.89400000000000e-05,2.97300000000000e-05,3.05900000000000e-05,3.15000000000000e-05,3.27100000000000e-05,3.39800000000000e-05,3.53000000000000e-05;2.62900000000000e-05,2.64600000000000e-05,2.65500000000000e-05,2.66400000000000e-05,2.67400000000000e-05,2.68400000000000e-05,2.69500000000000e-05,2.71800000000000e-05,2.74200000000000e-05,2.76300000000000e-05,2.76800000000000e-05,2.79600000000000e-05,2.85600000000000e-05,2.92300000000000e-05,2.99600000000000e-05,3.07500000000000e-05,3.15900000000000e-05,3.27100000000000e-05,3.38800000000000e-05,3.51100000000000e-05;2.68000000000000e-05,2.69600000000000e-05,2.70400000000000e-05,2.71300000000000e-05,2.72300000000000e-05,2.73200000000000e-05,2.74200000000000e-05,2.76300000000000e-05,2.78600000000000e-05,2.80500000000000e-05,2.81000000000000e-05,2.83600000000000e-05,2.89300000000000e-05,2.95500000000000e-05,3.02300000000000e-05,3.09600000000000e-05,3.17400000000000e-05,3.27700000000000e-05,3.38600000000000e-05,3.50000000000000e-05;2.73100000000000e-05,2.74600000000000e-05,2.75400000000000e-05,2.76200000000000e-05,2.77100000000000e-05,2.78000000000000e-05,2.78900000000000e-05,2.80900000000000e-05,2.83100000000000e-05,2.84900000000000e-05,2.85400000000000e-05,2.87800000000000e-05,2.93100000000000e-05,2.98900000000000e-05,3.05200000000000e-05,3.12000000000000e-05,3.19200000000000e-05,3.28900000000000e-05,3.39000000000000e-05,3.49600000000000e-05;2.78100000000000e-05,2.79600000000000e-05,2.80300000000000e-05,2.81100000000000e-05,2.81900000000000e-05,2.82800000000000e-05,2.83700000000000e-05,2.85600000000000e-05,2.87600000000000e-05,2.89300000000000e-05,2.89700000000000e-05,2.92000000000000e-05,2.96900000000000e-05,3.02400000000000e-05,3.08300000000000e-05,3.14700000000000e-05,3.21400000000000e-05,3.30400000000000e-05,3.39900000000000e-05,3.49900000000000e-05;2.83100000000000e-05,2.84500000000000e-05,2.85200000000000e-05,2.86000000000000e-05,2.86700000000000e-05,2.87500000000000e-05,2.88400000000000e-05,2.90200000000000e-05,2.92100000000000e-05,2.93700000000000e-05,2.94100000000000e-05,2.96300000000000e-05,3.00900000000000e-05,3.06000000000000e-05,3.11600000000000e-05,3.17500000000000e-05,3.23900000000000e-05,3.32300000000000e-05,3.41200000000000e-05,3.50600000000000e-05;2.88100000000000e-05,2.89400000000000e-05,2.90100000000000e-05,2.90800000000000e-05,2.91500000000000e-05,2.92300000000000e-05,2.93100000000000e-05,2.94800000000000e-05,2.96600000000000e-05,2.98100000000000e-05,2.98500000000000e-05,3.00500000000000e-05,3.04900000000000e-05,3.09800000000000e-05,3.15000000000000e-05,3.20600000000000e-05,3.26600000000000e-05,3.34500000000000e-05,3.42900000000000e-05,3.51700000000000e-05;2.93000000000000e-05,2.94300000000000e-05,2.94900000000000e-05,2.95600000000000e-05,2.96300000000000e-05,2.97000000000000e-05,2.97800000000000e-05,2.99400000000000e-05,3.01100000000000e-05,3.02600000000000e-05,3.02900000000000e-05,3.04900000000000e-05,3.09000000000000e-05,3.13600000000000e-05,3.18500000000000e-05,3.23800000000000e-05,3.29500000000000e-05,3.37000000000000e-05,3.44900000000000e-05,3.53200000000000e-05;2.97900000000000e-05,2.99100000000000e-05,2.99700000000000e-05,3.00400000000000e-05,3.01000000000000e-05,3.01700000000000e-05,3.02500000000000e-05,3.04000000000000e-05,3.05600000000000e-05,3.07000000000000e-05,3.07400000000000e-05,3.09200000000000e-05,3.13100000000000e-05,3.17500000000000e-05,3.22100000000000e-05,3.27200000000000e-05,3.32500000000000e-05,3.39600000000000e-05,3.47100000000000e-05,3.54900000000000e-05;3.02800000000000e-05,3.03900000000000e-05,3.04500000000000e-05,3.05100000000000e-05,3.05800000000000e-05,3.06400000000000e-05,3.07100000000000e-05,3.08600000000000e-05,3.10100000000000e-05,3.11400000000000e-05,3.11800000000000e-05,3.13500000000000e-05,3.17300000000000e-05,3.21400000000000e-05,3.25800000000000e-05,3.30600000000000e-05,3.35700000000000e-05,3.42400000000000e-05,3.49500000000000e-05,3.56900000000000e-05;3.07600000000000e-05,3.08600000000000e-05,3.09200000000000e-05,3.09800000000000e-05,3.10400000000000e-05,3.11100000000000e-05,3.11700000000000e-05,3.13100000000000e-05,3.14600000000000e-05,3.15900000000000e-05,3.16200000000000e-05,3.17900000000000e-05,3.21400000000000e-05,3.25400000000000e-05,3.29600000000000e-05,3.34100000000000e-05,3.38900000000000e-05,3.45300000000000e-05,3.52000000000000e-05,3.59100000000000e-05;3.12300000000000e-05,3.13400000000000e-05,3.13900000000000e-05,3.14500000000000e-05,3.15100000000000e-05,3.15700000000000e-05,3.16300000000000e-05,3.17700000000000e-05,3.19100000000000e-05,3.20300000000000e-05,3.20600000000000e-05,3.22200000000000e-05,3.25600000000000e-05,3.29300000000000e-05,3.33400000000000e-05,3.37700000000000e-05,3.42300000000000e-05,3.48300000000000e-05,3.54800000000000e-05,3.61500000000000e-05;3.17100000000000e-05,3.18100000000000e-05,3.18600000000000e-05,3.19100000000000e-05,3.19700000000000e-05,3.20300000000000e-05,3.20900000000000e-05,3.22200000000000e-05,3.23600000000000e-05,3.24700000000000e-05,3.25000000000000e-05,3.26500000000000e-05,3.29800000000000e-05,3.33400000000000e-05,3.37200000000000e-05,3.41300000000000e-05,3.45700000000000e-05,3.51500000000000e-05,3.57600000000000e-05,3.64100000000000e-05;3.21700000000000e-05,3.22700000000000e-05,3.23200000000000e-05,3.23700000000000e-05,3.24300000000000e-05,3.24900000000000e-05,3.25400000000000e-05,3.26700000000000e-05,3.28000000000000e-05,3.29100000000000e-05,3.29400000000000e-05,3.30800000000000e-05,3.34000000000000e-05,3.37400000000000e-05,3.41100000000000e-05,3.45000000000000e-05,3.49200000000000e-05,3.54700000000000e-05,3.60600000000000e-05,3.66700000000000e-05;3.26400000000000e-05,3.27300000000000e-05,3.27800000000000e-05,3.28300000000000e-05,3.28800000000000e-05,3.29400000000000e-05,3.30000000000000e-05,3.31100000000000e-05,3.32400000000000e-05,3.33500000000000e-05,3.33700000000000e-05,3.35100000000000e-05,3.38100000000000e-05,3.41400000000000e-05,3.44900000000000e-05,3.48700000000000e-05,3.52700000000000e-05,3.58000000000000e-05,3.63600000000000e-05,3.69500000000000e-05;3.31000000000000e-05,3.31900000000000e-05,3.32400000000000e-05,3.32800000000000e-05,3.33400000000000e-05,3.33900000000000e-05,3.34400000000000e-05,3.35600000000000e-05,3.36800000000000e-05,3.37800000000000e-05,3.38100000000000e-05,3.39400000000000e-05,3.42300000000000e-05,3.45400000000000e-05,3.48800000000000e-05,3.52400000000000e-05,3.56300000000000e-05,3.61300000000000e-05,3.66700000000000e-05,3.72400000000000e-05;3.35600000000000e-05,3.36400000000000e-05,3.36900000000000e-05,3.37400000000000e-05,3.37800000000000e-05,3.38400000000000e-05,3.38900000000000e-05,3.40000000000000e-05,3.41100000000000e-05,3.42100000000000e-05,3.42400000000000e-05,3.43700000000000e-05,3.46400000000000e-05,3.49500000000000e-05,3.52700000000000e-05,3.56200000000000e-05,3.59900000000000e-05,3.64700000000000e-05,3.69900000000000e-05,3.75300000000000e-05;3.40100000000000e-05,3.40900000000000e-05,3.41400000000000e-05,3.41800000000000e-05,3.42300000000000e-05,3.42800000000000e-05,3.43300000000000e-05,3.44400000000000e-05,3.45500000000000e-05,3.46400000000000e-05,3.46700000000000e-05,3.47900000000000e-05,3.50600000000000e-05,3.53500000000000e-05,3.56600000000000e-05,3.59900000000000e-05,3.63500000000000e-05,3.68200000000000e-05,3.73100000000000e-05,3.78300000000000e-05;3.44600000000000e-05,3.45400000000000e-05,3.45800000000000e-05,3.46200000000000e-05,3.46700000000000e-05,3.47200000000000e-05,3.47700000000000e-05,3.48700000000000e-05,3.49800000000000e-05,3.50700000000000e-05,3.50900000000000e-05,3.52100000000000e-05,3.54700000000000e-05,3.57500000000000e-05,3.60500000000000e-05,3.63700000000000e-05,3.67100000000000e-05,3.71600000000000e-05,3.76400000000000e-05,3.81400000000000e-05;3.49000000000000e-05,3.49800000000000e-05,3.50200000000000e-05,3.50600000000000e-05,3.51100000000000e-05,3.51500000000000e-05,3.52000000000000e-05,3.53000000000000e-05,3.54100000000000e-05,3.54900000000000e-05,3.55200000000000e-05,3.56300000000000e-05,3.58800000000000e-05,3.61500000000000e-05,3.64400000000000e-05,3.67500000000000e-05,3.70800000000000e-05,3.75100000000000e-05,3.79700000000000e-05,3.84500000000000e-05;3.53400000000000e-05,3.54200000000000e-05,3.54600000000000e-05,3.55000000000000e-05,3.55400000000000e-05,3.55900000000000e-05,3.56300000000000e-05,3.57300000000000e-05,3.58300000000000e-05,3.59200000000000e-05,3.59400000000000e-05,3.60500000000000e-05,3.62900000000000e-05,3.65500000000000e-05,3.68300000000000e-05,3.71300000000000e-05,3.74400000000000e-05,3.78600000000000e-05,3.83100000000000e-05,3.87700000000000e-05;3.57800000000000e-05,3.58500000000000e-05,3.58900000000000e-05,3.59300000000000e-05,3.59700000000000e-05,3.60200000000000e-05,3.60600000000000e-05,3.61500000000000e-05,3.62500000000000e-05,3.63400000000000e-05,3.63600000000000e-05,3.64600000000000e-05,3.67000000000000e-05,3.69500000000000e-05,3.72200000000000e-05,3.75100000000000e-05,3.78100000000000e-05,3.82200000000000e-05,3.86400000000000e-05,3.90900000000000e-05;3.62100000000000e-05,3.62800000000000e-05,3.63200000000000e-05,3.63600000000000e-05,3.64000000000000e-05,3.64400000000000e-05,3.64900000000000e-05,3.65800000000000e-05,3.66700000000000e-05,3.67500000000000e-05,3.67700000000000e-05,3.68800000000000e-05,3.71000000000000e-05,3.73500000000000e-05,3.76100000000000e-05,3.78900000000000e-05,3.81800000000000e-05,3.85700000000000e-05,3.89800000000000e-05,3.94200000000000e-05;3.66400000000000e-05,3.67100000000000e-05,3.67500000000000e-05,3.67900000000000e-05,3.68300000000000e-05,3.68700000000000e-05,3.69100000000000e-05,3.70000000000000e-05,3.70900000000000e-05,3.71700000000000e-05,3.71900000000000e-05,3.72900000000000e-05,3.75000000000000e-05,3.77400000000000e-05,3.79900000000000e-05,3.82600000000000e-05,3.85500000000000e-05,3.89300000000000e-05,3.93200000000000e-05,3.97400000000000e-05;3.70700000000000e-05,3.71400000000000e-05,3.71700000000000e-05,3.72100000000000e-05,3.72500000000000e-05,3.72900000000000e-05,3.73300000000000e-05,3.74100000000000e-05,3.75000000000000e-05,3.75800000000000e-05,3.76000000000000e-05,3.76900000000000e-05,3.79100000000000e-05,3.81300000000000e-05,3.83800000000000e-05,3.86400000000000e-05,3.89200000000000e-05,3.92800000000000e-05,3.96700000000000e-05,4.00700000000000e-05;3.74900000000000e-05,3.75600000000000e-05,3.75900000000000e-05,3.76300000000000e-05,3.76700000000000e-05,3.77000000000000e-05,3.77400000000000e-05,3.78300000000000e-05,3.79100000000000e-05,3.79900000000000e-05,3.80000000000000e-05,3.81000000000000e-05,3.83000000000000e-05,3.85300000000000e-05,3.87600000000000e-05,3.90200000000000e-05,3.92800000000000e-05,3.96400000000000e-05,4.00100000000000e-05,4.04000000000000e-05;3.79100000000000e-05,3.79800000000000e-05,3.80100000000000e-05,3.80400000000000e-05,3.80800000000000e-05,3.81200000000000e-05,3.81600000000000e-05,3.82400000000000e-05,3.83200000000000e-05,3.83900000000000e-05,3.84100000000000e-05,3.85000000000000e-05,3.87000000000000e-05,3.89100000000000e-05,3.91500000000000e-05,3.93900000000000e-05,3.96500000000000e-05,3.99900000000000e-05,4.03600000000000e-05,4.07400000000000e-05;3.83300000000000e-05,3.83900000000000e-05,3.84200000000000e-05,3.84600000000000e-05,3.84900000000000e-05,3.85300000000000e-05,3.85700000000000e-05,3.86400000000000e-05,3.87200000000000e-05,3.87900000000000e-05,3.88100000000000e-05,3.89000000000000e-05,3.90900000000000e-05,3.93000000000000e-05,3.95300000000000e-05,3.97600000000000e-05,4.00200000000000e-05,4.03500000000000e-05,4.07000000000000e-05,4.10700000000000e-05;3.87400000000000e-05,3.88000000000000e-05,3.88300000000000e-05,3.88700000000000e-05,3.89000000000000e-05,3.89400000000000e-05,3.89700000000000e-05,3.90500000000000e-05,3.91300000000000e-05,3.91900000000000e-05,3.92100000000000e-05,3.93000000000000e-05,3.94900000000000e-05,3.96900000000000e-05,3.99100000000000e-05,4.01400000000000e-05,4.03800000000000e-05,4.07000000000000e-05,4.10400000000000e-05,4.14000000000000e-05;3.91500000000000e-05,3.92100000000000e-05,3.92400000000000e-05,3.92700000000000e-05,3.93100000000000e-05,3.93400000000000e-05,3.93700000000000e-05,3.94500000000000e-05,3.95300000000000e-05,3.95900000000000e-05,3.96100000000000e-05,3.96900000000000e-05,3.98800000000000e-05,4.00700000000000e-05,4.02800000000000e-05,4.05100000000000e-05,4.07400000000000e-05,4.10600000000000e-05,4.13900000000000e-05,4.17400000000000e-05;3.95600000000000e-05,3.96100000000000e-05,3.96400000000000e-05,3.96700000000000e-05,3.97100000000000e-05,3.97400000000000e-05,3.97800000000000e-05,3.98500000000000e-05,3.99200000000000e-05,3.99900000000000e-05,4.00000000000000e-05,4.00900000000000e-05,4.02600000000000e-05,4.04500000000000e-05,4.06600000000000e-05,4.08800000000000e-05,4.11100000000000e-05,4.14100000000000e-05,4.17300000000000e-05,4.20700000000000e-05;3.99600000000000e-05,4.00100000000000e-05,4.00400000000000e-05,4.00700000000000e-05,4.01100000000000e-05,4.01400000000000e-05,4.01700000000000e-05,4.02400000000000e-05,4.03200000000000e-05,4.03800000000000e-05,4.03900000000000e-05,4.04700000000000e-05,4.06500000000000e-05,4.08300000000000e-05,4.10300000000000e-05,4.12400000000000e-05,4.14700000000000e-05,4.17700000000000e-05,4.20800000000000e-05,4.24100000000000e-05;4.03600000000000e-05,4.04100000000000e-05,4.04400000000000e-05,4.04700000000000e-05,4.05000000000000e-05,4.05300000000000e-05,4.05700000000000e-05,4.06400000000000e-05,4.07100000000000e-05,4.07700000000000e-05,4.07800000000000e-05,4.08600000000000e-05,4.10300000000000e-05,4.12100000000000e-05,4.14100000000000e-05,4.16100000000000e-05,4.18300000000000e-05,4.21200000000000e-05,4.24200000000000e-05,4.27400000000000e-05;4.07500000000000e-05,4.08100000000000e-05,4.08400000000000e-05,4.08700000000000e-05,4.09000000000000e-05,4.09300000000000e-05,4.09600000000000e-05,4.10300000000000e-05,4.11000000000000e-05,4.11500000000000e-05,4.11700000000000e-05,4.12500000000000e-05,4.14100000000000e-05,4.15900000000000e-05,4.17800000000000e-05,4.19800000000000e-05,4.21900000000000e-05,4.24700000000000e-05,4.27700000000000e-05,4.30800000000000e-05];
    y = y';
    val = interpinbounds(T,P,y,temp,press);
    mu = val;
end

function k = CO2_k_real(pos,time,temp,press)
    %for 7.5 < P < 32.5 MPa, 20 < T < 700 C
    T = [20;33.6000000000000;47.2000000000000;60.8000000000000;74.4000000000000;88;101.600000000000;115.200000000000;128.800000000000;142.400000000000;156;169.600000000000;183.200000000000;196.800000000000;210.400000000000;224;237.600000000000;251.200000000000;264.800000000000;278.400000000000;292;305.600000000000;319.200000000000;332.800000000000;346.400000000000;360;373.600000000000;387.200000000000;400.800000000000;414.400000000000;428;441.600000000000;455.200000000000;468.800000000000;482.400000000000;496;509.600000000000;523.200000000000;536.800000000000;550.400000000000;564;577.600000000000;591.200000000000;604.800000000000;618.400000000000;632;645.600000000000;659.200000000000;672.800000000000;686.400000000000;700];
    P = [7.50000000000000,8.50000000000000,9,9.50000000000000,10,10.5000000000000,11,12,13,13.8000000000000,14,15,17,19,21,23,25,27.5000000000000,30,32.5000000000000];
    P = P*1e6;
    [T,P] = meshgrid(T,P);
    y = [0.0882700000000000,0.0914000000000000,0.0927900000000000,0.0940900000000000,0.0953100000000000,8.27600000000000e-05,0.0975900000000000,0.0996700000000000,0.101600000000000,0.103000000000000,0.103400000000000,0.105100000000000,0.108300000000000,0.111200000000000,0.113800000000000,0.116400000000000,0.118700000000000,0.121500000000000,0.124200000000000,0.126700000000000;0.0513200000000000,0.0750100000000000,0.0758600000000000,0.0773300000000000,0.0788100000000000,6.22600000000000e-05,0.0810500000000000,0.0825700000000000,0.0854800000000000,0.0875600000000000,0.0880500000000000,0.0903700000000000,0.0944600000000000,0.0980400000000000,0.101300000000000,0.104200000000000,0.106900000000000,0.110100000000000,0.113000000000000,0.115800000000000;0.0307000000000000,0.0383300000000000,0.0443900000000000,0.0521200000000000,0.0595100000000000,3.69200000000000e-05,0.0654600000000000,0.0686800000000000,0.0717900000000000,0.0739300000000000,0.0744000000000000,0.0750000000000000,0.0806300000000000,0.0852100000000000,0.0891300000000000,0.0926100000000000,0.0957600000000000,0.0993300000000000,0.102600000000000,0.105600000000000;0.0280600000000000,0.0313600000000000,0.0334700000000000,0.0359500000000000,0.0388200000000000,2.51600000000000e-05,0.0455400000000000,0.0522900000000000,0.0574600000000000,0.0605900000000000,0.0612900000000000,0.0644900000000000,0.0696300000000000,0.0728700000000000,0.0776500000000000,0.0817400000000000,0.0853600000000000,0.0893900000000000,0.0930100000000000,0.0963200000000000;0.0274400000000000,0.0295300000000000,0.0307600000000000,0.0321400000000000,0.0336800000000000,2.29900000000000e-05,0.0372100000000000,0.0412700000000000,0.0455800000000000,0.0489400000000000,0.0497400000000000,0.0534700000000000,0.0596600000000000,0.0645700000000000,0.0672700000000000,0.0718900000000000,0.0759300000000000,0.0803700000000000,0.0843200000000000,0.0878900000000000;0.0275400000000000,0.0290700000000000,0.0299400000000000,0.0308800000000000,0.0319100000000000,2.23400000000000e-05,0.0342200000000000,0.0368400000000000,0.0397300000000000,0.0421700000000000,0.0427900000000000,0.0458700000000000,0.0517000000000000,0.0568100000000000,0.0611400000000000,0.0641400000000000,0.0677900000000000,0.0724800000000000,0.0766500000000000,0.0804100000000000;0.0279800000000000,0.0291900000000000,0.0298600000000000,0.0305800000000000,0.0313600000000000,2.22100000000000e-05,0.0330500000000000,0.0349500000000000,0.0370300000000000,0.0388100000000000,0.0392700000000000,0.0416200000000000,0.0464100000000000,0.0510200000000000,0.0552500000000000,0.0589900000000000,0.0617700000000000,0.0659400000000000,0.0701400000000000,0.0739700000000000;0.0286100000000000,0.0296200000000000,0.0301700000000000,0.0307600000000000,0.0313700000000000,2.23100000000000e-05,0.0327100000000000,0.0341900000000000,0.0358000000000000,0.0371700000000000,0.0375300000000000,0.0393500000000000,0.0432000000000000,0.0471200000000000,0.0509200000000000,0.0544900000000000,0.0577100000000000,0.0607900000000000,0.0648700000000000,0.0686300000000000;0.0293600000000000,0.0302300000000000,0.0307000000000000,0.0311900000000000,0.0317100000000000,2.25400000000000e-05,0.0328100000000000,0.0340200000000000,0.0353300000000000,0.0364300000000000,0.0367200000000000,0.0382000000000000,0.0413400000000000,0.0446300000000000,0.0479400000000000,0.0511600000000000,0.0542200000000000,0.0576600000000000,0.0607500000000000,0.0643400000000000;0.0301900000000000,0.0309500000000000,0.0313600000000000,0.0317900000000000,0.0322300000000000,2.28600000000000e-05,0.0331800000000000,0.0342000000000000,0.0352900000000000,0.0362200000000000,0.0364600000000000,0.0376800000000000,0.0403100000000000,0.0431000000000000,0.0459600000000000,0.0488200000000000,0.0516200000000000,0.0549400000000000,0.0576300000000000,0.0609900000000000;0.0310700000000000,0.0317500000000000,0.0321100000000000,0.0324900000000000,0.0328800000000000,2.32200000000000e-05,0.0337000000000000,0.0345800000000000,0.0355200000000000,0.0363100000000000,0.0365200000000000,0.0375700000000000,0.0398000000000000,0.0421900000000000,0.0446800000000000,0.0472100000000000,0.0497300000000000,0.0528200000000000,0.0557300000000000,0.0584400000000000;0.0319800000000000,0.0325900000000000,0.0329200000000000,0.0332500000000000,0.0335900000000000,2.36300000000000e-05,0.0343200000000000,0.0350900000000000,0.0359100000000000,0.0366000000000000,0.0367700000000000,0.0376800000000000,0.0396100000000000,0.0416800000000000,0.0438500000000000,0.0460900000000000,0.0483600000000000,0.0511900000000000,0.0539700000000000,0.0565500000000000;0.0329100000000000,0.0334600000000000,0.0337500000000000,0.0340500000000000,0.0343600000000000,2.40500000000000e-05,0.0350000000000000,0.0356800000000000,0.0364000000000000,0.0370000000000000,0.0371500000000000,0.0379400000000000,0.0396100000000000,0.0414100000000000,0.0433100000000000,0.0453000000000000,0.0473400000000000,0.0499400000000000,0.0525600000000000,0.0551700000000000;0.0339100000000000,0.0344200000000000,0.0346900000000000,0.0349700000000000,0.0352500000000000,2.45000000000000e-05,0.0358500000000000,0.0364800000000000,0.0371400000000000,0.0376800000000000,0.0378200000000000,0.0385400000000000,0.0400700000000000,0.0417000000000000,0.0434100000000000,0.0452000000000000,0.0470500000000000,0.0494100000000000,0.0518100000000000,0.0542000000000000;0.0349100000000000,0.0353900000000000,0.0356500000000000,0.0359100000000000,0.0361700000000000,2.49500000000000e-05,0.0367300000000000,0.0373100000000000,0.0379200000000000,0.0384300000000000,0.0385600000000000,0.0392200000000000,0.0406200000000000,0.0421100000000000,0.0436800000000000,0.0453100000000000,0.0470000000000000,0.0491600000000000,0.0513500000000000,0.0535600000000000;0.0359200000000000,0.0363800000000000,0.0366100000000000,0.0368600000000000,0.0371100000000000,2.54200000000000e-05,0.0376300000000000,0.0381700000000000,0.0387400000000000,0.0392200000000000,0.0393400000000000,0.0399500000000000,0.0412400000000000,0.0426200000000000,0.0440600000000000,0.0455600000000000,0.0471100000000000,0.0491000000000000,0.0511300000000000,0.0531700000000000;0.0369300000000000,0.0373700000000000,0.0375900000000000,0.0378200000000000,0.0380600000000000,2.58900000000000e-05,0.0385500000000000,0.0390600000000000,0.0395900000000000,0.0400400000000000,0.0401500000000000,0.0407200000000000,0.0419300000000000,0.0432000000000000,0.0445400000000000,0.0459300000000000,0.0473600000000000,0.0492100000000000,0.0510900000000000,0.0529900000000000;0.0379500000000000,0.0383600000000000,0.0385700000000000,0.0387900000000000,0.0390200000000000,2.63700000000000e-05,0.0394800000000000,0.0399600000000000,0.0404700000000000,0.0408800000000000,0.0409900000000000,0.0415300000000000,0.0426500000000000,0.0438400000000000,0.0450900000000000,0.0463900000000000,0.0477200000000000,0.0494400000000000,0.0511900000000000,0.0529700000000000;0.0389600000000000,0.0393600000000000,0.0395600000000000,0.0397700000000000,0.0399800000000000,2.68400000000000e-05,0.0404200000000000,0.0408800000000000,0.0413600000000000,0.0417500000000000,0.0418500000000000,0.0423500000000000,0.0434200000000000,0.0445300000000000,0.0457000000000000,0.0469100000000000,0.0481600000000000,0.0497700000000000,0.0514100000000000,0.0530800000000000;0.0399800000000000,0.0403600000000000,0.0405500000000000,0.0407500000000000,0.0409600000000000,2.73200000000000e-05,0.0413700000000000,0.0418100000000000,0.0422600000000000,0.0426300000000000,0.0427200000000000,0.0432000000000000,0.0442100000000000,0.0452600000000000,0.0463600000000000,0.0475000000000000,0.0486800000000000,0.0501800000000000,0.0517300000000000,0.0533000000000000;0.0410000000000000,0.0413600000000000,0.0415400000000000,0.0417400000000000,0.0419300000000000,2.78000000000000e-05,0.0423300000000000,0.0427400000000000,0.0431700000000000,0.0435300000000000,0.0436200000000000,0.0440700000000000,0.0450200000000000,0.0460200000000000,0.0470600000000000,0.0481300000000000,0.0492400000000000,0.0506700000000000,0.0521200000000000,0.0536000000000000;0.0420100000000000,0.0423600000000000,0.0425400000000000,0.0427200000000000,0.0429100000000000,2.82800000000000e-05,0.0432900000000000,0.0436900000000000,0.0440900000000000,0.0444300000000000,0.0445200000000000,0.0449500000000000,0.0458500000000000,0.0468000000000000,0.0477800000000000,0.0488000000000000,0.0498500000000000,0.0512000000000000,0.0525800000000000,0.0539800000000000;0.0430200000000000,0.0433600000000000,0.0435300000000000,0.0437100000000000,0.0438800000000000,2.87500000000000e-05,0.0442500000000000,0.0446300000000000,0.0450200000000000,0.0453400000000000,0.0454300000000000,0.0458400000000000,0.0467000000000000,0.0476000000000000,0.0485400000000000,0.0495100000000000,0.0505000000000000,0.0517800000000000,0.0530900000000000,0.0544300000000000;0.0440300000000000,0.0443600000000000,0.0445200000000000,0.0446900000000000,0.0448600000000000,2.92300000000000e-05,0.0452200000000000,0.0455800000000000,0.0459500000000000,0.0462600000000000,0.0463400000000000,0.0467400000000000,0.0475600000000000,0.0484200000000000,0.0493100000000000,0.0502400000000000,0.0511900000000000,0.0524000000000000,0.0536500000000000,0.0549200000000000;0.0450400000000000,0.0453500000000000,0.0455100000000000,0.0456800000000000,0.0458400000000000,2.97000000000000e-05,0.0461800000000000,0.0465300000000000,0.0468900000000000,0.0471900000000000,0.0472600000000000,0.0476400000000000,0.0484300000000000,0.0492500000000000,0.0501100000000000,0.0509900000000000,0.0518900000000000,0.0530600000000000,0.0542500000000000,0.0554600000000000;0.0460500000000000,0.0463500000000000,0.0465000000000000,0.0466600000000000,0.0468200000000000,3.01700000000000e-05,0.0471400000000000,0.0474800000000000,0.0478300000000000,0.0481100000000000,0.0481800000000000,0.0485500000000000,0.0493100000000000,0.0500900000000000,0.0509100000000000,0.0517600000000000,0.0526200000000000,0.0537400000000000,0.0548700000000000,0.0560300000000000;0.0470500000000000,0.0473400000000000,0.0474900000000000,0.0476400000000000,0.0477900000000000,3.06400000000000e-05,0.0481100000000000,0.0484300000000000,0.0487700000000000,0.0490400000000000,0.0491100000000000,0.0494600000000000,0.0501900000000000,0.0509500000000000,0.0517300000000000,0.0525400000000000,0.0533700000000000,0.0544400000000000,0.0555300000000000,0.0566400000000000;0.0480400000000000,0.0483300000000000,0.0484700000000000,0.0486200000000000,0.0487700000000000,3.11100000000000e-05,0.0490700000000000,0.0493800000000000,0.0497100000000000,0.0499700000000000,0.0500400000000000,0.0503800000000000,0.0510800000000000,0.0518100000000000,0.0525600000000000,0.0533400000000000,0.0541400000000000,0.0551600000000000,0.0562100000000000,0.0572700000000000;0.0490400000000000,0.0493100000000000,0.0494500000000000,0.0495900000000000,0.0497400000000000,3.15700000000000e-05,0.0500300000000000,0.0503400000000000,0.0506500000000000,0.0509000000000000,0.0509700000000000,0.0512900000000000,0.0519700000000000,0.0526700000000000,0.0534000000000000,0.0541500000000000,0.0549100000000000,0.0559000000000000,0.0569100000000000,0.0579300000000000;0.0500300000000000,0.0502900000000000,0.0504300000000000,0.0505600000000000,0.0507000000000000,3.20300000000000e-05,0.0509900000000000,0.0512800000000000,0.0515900000000000,0.0518300000000000,0.0519000000000000,0.0522100000000000,0.0528600000000000,0.0535400000000000,0.0542400000000000,0.0549600000000000,0.0557000000000000,0.0566500000000000,0.0576200000000000,0.0586100000000000;0.0510100000000000,0.0512700000000000,0.0514000000000000,0.0515300000000000,0.0516700000000000,3.24900000000000e-05,0.0519500000000000,0.0522300000000000,0.0525200000000000,0.0527600000000000,0.0528200000000000,0.0531300000000000,0.0537600000000000,0.0544200000000000,0.0550900000000000,0.0557900000000000,0.0565000000000000,0.0574200000000000,0.0583500000000000,0.0593100000000000;0.0519900000000000,0.0522400000000000,0.0523700000000000,0.0525000000000000,0.0526300000000000,3.29400000000000e-05,0.0529000000000000,0.0531800000000000,0.0534600000000000,0.0536900000000000,0.0537500000000000,0.0540500000000000,0.0546600000000000,0.0552900000000000,0.0559500000000000,0.0566200000000000,0.0573100000000000,0.0581900000000000,0.0591000000000000,0.0600200000000000;0.0529700000000000,0.0532100000000000,0.0533400000000000,0.0534600000000000,0.0535900000000000,3.33900000000000e-05,0.0538500000000000,0.0541200000000000,0.0544000000000000,0.0546200000000000,0.0546800000000000,0.0549700000000000,0.0555600000000000,0.0561700000000000,0.0568100000000000,0.0574600000000000,0.0581200000000000,0.0589800000000000,0.0598500000000000,0.0607400000000000;0.0539400000000000,0.0541800000000000,0.0543000000000000,0.0544200000000000,0.0545500000000000,3.38400000000000e-05,0.0548000000000000,0.0550600000000000,0.0553300000000000,0.0555500000000000,0.0556000000000000,0.0558800000000000,0.0564600000000000,0.0570500000000000,0.0576700000000000,0.0583000000000000,0.0589400000000000,0.0597700000000000,0.0606200000000000,0.0614800000000000;0.0549100000000000,0.0551400000000000,0.0552600000000000,0.0553800000000000,0.0555000000000000,3.42800000000000e-05,0.0557500000000000,0.0560000000000000,0.0562600000000000,0.0564700000000000,0.0565300000000000,0.0568000000000000,0.0573600000000000,0.0579300000000000,0.0585300000000000,0.0591400000000000,0.0597700000000000,0.0605700000000000,0.0613900000000000,0.0622200000000000;0.0558700000000000,0.0561000000000000,0.0562100000000000,0.0563300000000000,0.0564500000000000,3.47200000000000e-05,0.0566900000000000,0.0569400000000000,0.0571900000000000,0.0574000000000000,0.0574500000000000,0.0577100000000000,0.0582500000000000,0.0588200000000000,0.0593900000000000,0.0599900000000000,0.0606000000000000,0.0613700000000000,0.0621700000000000,0.0629800000000000;0.0568300000000000,0.0570500000000000,0.0571600000000000,0.0572800000000000,0.0573900000000000,3.51500000000000e-05,0.0576300000000000,0.0578700000000000,0.0581200000000000,0.0583200000000000,0.0583700000000000,0.0586200000000000,0.0591500000000000,0.0597000000000000,0.0602600000000000,0.0608400000000000,0.0614300000000000,0.0621800000000000,0.0629600000000000,0.0637400000000000;0.0577800000000000,0.0580000000000000,0.0581100000000000,0.0582200000000000,0.0583300000000000,3.55900000000000e-05,0.0585600000000000,0.0588000000000000,0.0590400000000000,0.0592300000000000,0.0592800000000000,0.0595300000000000,0.0600500000000000,0.0605800000000000,0.0611300000000000,0.0616900000000000,0.0622600000000000,0.0630000000000000,0.0637500000000000,0.0645100000000000;0.0587300000000000,0.0589400000000000,0.0590500000000000,0.0591600000000000,0.0592700000000000,3.60200000000000e-05,0.0595000000000000,0.0597200000000000,0.0599600000000000,0.0601500000000000,0.0602000000000000,0.0604400000000000,0.0609400000000000,0.0614600000000000,0.0619900000000000,0.0625400000000000,0.0631000000000000,0.0638100000000000,0.0645400000000000,0.0652800000000000;0.0596800000000000,0.0598800000000000,0.0599900000000000,0.0601000000000000,0.0602000000000000,3.64400000000000e-05,0.0604200000000000,0.0606500000000000,0.0608800000000000,0.0610600000000000,0.0611100000000000,0.0613500000000000,0.0618400000000000,0.0623400000000000,0.0628600000000000,0.0633900000000000,0.0639300000000000,0.0646300000000000,0.0653400000000000,0.0660600000000000;0.0606200000000000,0.0608200000000000,0.0609200000000000,0.0610300000000000,0.0611300000000000,3.68700000000000e-05,0.0613500000000000,0.0615700000000000,0.0617900000000000,0.0619700000000000,0.0620200000000000,0.0622500000000000,0.0627300000000000,0.0632200000000000,0.0637200000000000,0.0642400000000000,0.0647700000000000,0.0654500000000000,0.0661400000000000,0.0668500000000000;0.0615600000000000,0.0617500000000000,0.0618600000000000,0.0619600000000000,0.0620600000000000,3.72900000000000e-05,0.0622700000000000,0.0624800000000000,0.0627000000000000,0.0628800000000000,0.0629200000000000,0.0631500000000000,0.0636200000000000,0.0641000000000000,0.0645900000000000,0.0650900000000000,0.0656100000000000,0.0662700000000000,0.0669500000000000,0.0676300000000000;0.0624900000000000,0.0626800000000000,0.0627800000000000,0.0628800000000000,0.0629800000000000,3.77000000000000e-05,0.0631900000000000,0.0634000000000000,0.0636100000000000,0.0637800000000000,0.0638300000000000,0.0640500000000000,0.0645000000000000,0.0649700000000000,0.0654500000000000,0.0659500000000000,0.0664500000000000,0.0670900000000000,0.0677500000000000,0.0684200000000000;0.0634200000000000,0.0636100000000000,0.0637000000000000,0.0638000000000000,0.0639000000000000,3.81200000000000e-05,0.0641000000000000,0.0643100000000000,0.0645100000000000,0.0646800000000000,0.0647300000000000,0.0649400000000000,0.0653900000000000,0.0658500000000000,0.0663200000000000,0.0668000000000000,0.0672900000000000,0.0679200000000000,0.0685600000000000,0.0692100000000000;0.0643400000000000,0.0645300000000000,0.0646200000000000,0.0647200000000000,0.0648100000000000,3.85300000000000e-05,0.0650100000000000,0.0652100000000000,0.0654200000000000,0.0655800000000000,0.0656200000000000,0.0658400000000000,0.0662700000000000,0.0667200000000000,0.0671800000000000,0.0676500000000000,0.0681300000000000,0.0687400000000000,0.0693700000000000,0.0700100000000000;0.0652600000000000,0.0654400000000000,0.0655400000000000,0.0656300000000000,0.0657300000000000,3.89400000000000e-05,0.0659200000000000,0.0661100000000000,0.0663100000000000,0.0664800000000000,0.0665200000000000,0.0667300000000000,0.0671500000000000,0.0675900000000000,0.0680400000000000,0.0685000000000000,0.0689700000000000,0.0695700000000000,0.0701800000000000,0.0708000000000000;0.0661800000000000,0.0663600000000000,0.0664500000000000,0.0665400000000000,0.0666300000000000,3.93400000000000e-05,0.0668200000000000,0.0670100000000000,0.0672100000000000,0.0673700000000000,0.0674100000000000,0.0676100000000000,0.0680300000000000,0.0684600000000000,0.0689000000000000,0.0693500000000000,0.0698100000000000,0.0703900000000000,0.0709900000000000,0.0716000000000000;0.0670900000000000,0.0672600000000000,0.0673500000000000,0.0674400000000000,0.0675400000000000,3.97400000000000e-05,0.0677200000000000,0.0679100000000000,0.0681000000000000,0.0682600000000000,0.0683000000000000,0.0685000000000000,0.0689000000000000,0.0693200000000000,0.0697500000000000,0.0701900000000000,0.0706400000000000,0.0712200000000000,0.0718000000000000,0.0724000000000000;0.0680000000000000,0.0681700000000000,0.0682600000000000,0.0683400000000000,0.0684300000000000,4.01400000000000e-05,0.0686200000000000,0.0688000000000000,0.0689900000000000,0.0691400000000000,0.0691800000000000,0.0693800000000000,0.0697800000000000,0.0701900000000000,0.0706100000000000,0.0710400000000000,0.0714800000000000,0.0720400000000000,0.0726100000000000,0.0732000000000000;0.0689000000000000,0.0690700000000000,0.0691500000000000,0.0692400000000000,0.0693300000000000,4.05300000000000e-05,0.0695100000000000,0.0696900000000000,0.0698700000000000,0.0700200000000000,0.0700600000000000,0.0702500000000000,0.0706400000000000,0.0710500000000000,0.0714600000000000,0.0718800000000000,0.0723100000000000,0.0728600000000000,0.0734200000000000,0.0740000000000000;0.0698000000000000,0.0699700000000000,0.0700500000000000,0.0701300000000000,0.0702200000000000,4.09300000000000e-05,0.0704000000000000,0.0705700000000000,0.0707600000000000,0.0709000000000000,0.0709400000000000,0.0711300000000000,0.0715100000000000,0.0719100000000000,0.0723100000000000,0.0727200000000000,0.0731500000000000,0.0736800000000000,0.0742300000000000,0.0747900000000000];
    y = y';
    val = interpinbounds(T,P,y,temp,press);
    k = val;
end

function cp = CO2_cp_real(pos,time,temp,press)
    %for 7.5 < P < 32.5 MPa, 20 < T < 700 C
    T = [20;33.6000000000000;47.2000000000000;60.8000000000000;74.4000000000000;88;101.600000000000;115.200000000000;128.800000000000;142.400000000000;156;169.600000000000;183.200000000000;196.800000000000;210.400000000000;224;237.600000000000;251.200000000000;264.800000000000;278.400000000000;292;305.600000000000;319.200000000000;332.800000000000;346.400000000000;360;373.600000000000;387.200000000000;400.800000000000;414.400000000000;428;441.600000000000;455.200000000000;468.800000000000;482.400000000000;496;509.600000000000;523.200000000000;536.800000000000;550.400000000000;564;577.600000000000;591.200000000000;604.800000000000;618.400000000000;632;645.600000000000;659.200000000000;672.800000000000;686.400000000000;700];
    P = [7.50000000000000,8.50000000000000,9,9.50000000000000,10,10.5000000000000,11,12,13,13.8000000000000,14,15,17,19,21,23,25,27.5000000000000,30,32.5000000000000];
    P = P*1e6;
    [T,P] = meshgrid(T,P);
    y = [3116,2861,2768,2689,2622,2563,2512,2425,2355,2307,2296,2246,2165,2102,2051,2008,1972,1934,1901,1874;8942,6536,4889,4155,3727,3440,3231,2943,2749,2633,2608,2499,2342,2232,2149,2085,2033,1980,1937,1901;2359,3541,4665,6252,7442,6901,5789,4320,3609,3267,3199,2933,2605,2409,2276,2179,2105,2033,1977,1931;1744,2089,2319,2599,2934,3315,3708,4270,4220,3918,3836,3454,2928,2621,2422,2284,2182,2087,2015,1959;1499,1676,1782,1900,2032,2177,2334,2666,2967,3130,3156,3196,2997,2733,2522,2365,2246,2133,2048,1982;1368,1481,1544,1612,1686,1764,1848,2026,2209,2348,2380,2517,2639,2591,2476,2357,2253,2147,2062,1993;1288,1368,1411,1457,1505,1556,1609,1720,1837,1930,1953,2062,2229,2305,2301,2252,2189,2110,2040,1980;1235,1295,1328,1361,1396,1432,1469,1547,1627,1693,1709,1789,1931,2032,2084,2094,2075,2034,1986,1940;1199,1247,1272,1297,1324,1351,1379,1437,1497,1545,1557,1617,1730,1824,1890,1929,1942,1935,1912,1883;1174,1213,1233,1253,1274,1296,1318,1363,1410,1447,1456,1503,1593,1673,1738,1784,1814,1829,1827,1815;1156,1188,1205,1222,1239,1257,1275,1311,1348,1378,1386,1423,1496,1563,1622,1668,1703,1730,1743,1744;1143,1171,1185,1199,1214,1228,1243,1273,1304,1329,1335,1366,1426,1483,1534,1577,1612,1644,1665,1675;1134,1158,1170,1182,1195,1207,1220,1245,1271,1292,1297,1323,1374,1422,1467,1506,1539,1572,1596,1612;1128,1149,1160,1170,1181,1192,1203,1225,1247,1265,1269,1291,1335,1376,1415,1450,1481,1513,1538,1557;1125,1143,1152,1161,1171,1180,1190,1209,1228,1244,1247,1267,1304,1341,1375,1406,1434,1464,1489,1509;1123,1139,1147,1155,1164,1172,1180,1197,1214,1228,1231,1248,1281,1313,1343,1371,1397,1425,1448,1468;1122,1136,1144,1151,1159,1166,1174,1189,1204,1216,1219,1233,1263,1291,1318,1343,1366,1392,1415,1434;1122,1135,1142,1149,1155,1162,1169,1182,1196,1206,1209,1222,1248,1274,1298,1321,1342,1366,1387,1405;1123,1135,1141,1147,1153,1159,1166,1178,1190,1199,1202,1214,1237,1260,1282,1303,1322,1344,1364,1381;1125,1136,1142,1147,1153,1158,1164,1175,1186,1194,1197,1207,1229,1249,1269,1288,1306,1326,1344,1361;1128,1138,1143,1148,1153,1158,1163,1173,1183,1191,1193,1203,1222,1241,1259,1276,1292,1311,1328,1344;1131,1140,1144,1149,1154,1158,1163,1172,1181,1189,1190,1199,1217,1234,1251,1267,1282,1299,1315,1330;1134,1142,1147,1151,1155,1159,1164,1172,1180,1187,1189,1197,1213,1229,1245,1259,1273,1289,1304,1318;1137,1145,1149,1153,1157,1161,1165,1173,1180,1187,1188,1196,1211,1226,1240,1253,1266,1281,1295,1308;1141,1148,1152,1156,1159,1163,1167,1174,1181,1187,1188,1195,1209,1223,1236,1248,1260,1275,1288,1300;1145,1152,1155,1158,1162,1165,1169,1175,1182,1188,1189,1195,1208,1221,1233,1245,1256,1269,1282,1293;1149,1155,1158,1162,1165,1168,1171,1177,1184,1189,1190,1196,1208,1220,1231,1242,1253,1265,1277,1288;1153,1159,1162,1165,1168,1171,1174,1180,1186,1190,1191,1197,1208,1219,1230,1240,1250,1262,1273,1283;1157,1163,1165,1168,1171,1174,1177,1182,1188,1192,1193,1199,1209,1220,1229,1239,1248,1259,1270,1279;1161,1167,1169,1172,1174,1177,1180,1185,1190,1194,1195,1200,1210,1220,1229,1238,1247,1258,1267,1276;1165,1170,1173,1175,1178,1180,1183,1188,1193,1197,1198,1202,1212,1221,1230,1238,1246,1256,1266,1274;1170,1174,1177,1179,1182,1184,1186,1191,1195,1199,1200,1205,1213,1222,1230,1238,1246,1256,1264,1273;1174,1178,1181,1183,1185,1187,1190,1194,1198,1202,1203,1207,1215,1223,1231,1239,1246,1255,1264,1271;1178,1182,1185,1187,1189,1191,1193,1197,1201,1205,1205,1209,1217,1225,1233,1240,1247,1255,1263,1271;1182,1186,1188,1191,1193,1195,1197,1200,1204,1208,1208,1212,1220,1227,1234,1241,1248,1256,1263,1270;1187,1190,1192,1194,1196,1198,1200,1204,1208,1210,1211,1215,1222,1229,1236,1242,1249,1256,1263,1270;1191,1194,1196,1198,1200,1202,1204,1207,1211,1213,1214,1218,1224,1231,1238,1244,1250,1257,1264,1271;1195,1198,1200,1202,1204,1205,1207,1210,1214,1217,1217,1221,1227,1233,1239,1245,1251,1258,1265,1271;1199,1202,1204,1206,1207,1209,1211,1214,1217,1220,1220,1223,1230,1236,1242,1247,1253,1259,1266,1272;1203,1206,1208,1209,1211,1213,1214,1217,1220,1223,1223,1226,1232,1238,1244,1249,1254,1261,1267,1273;1207,1210,1212,1213,1215,1216,1218,1221,1224,1226,1226,1229,1235,1241,1246,1251,1256,1262,1268,1274;1211,1214,1215,1217,1218,1220,1221,1224,1227,1229,1230,1232,1238,1243,1248,1253,1258,1264,1269,1275;1215,1217,1219,1220,1222,1223,1224,1227,1230,1232,1233,1235,1240,1246,1250,1255,1260,1266,1271,1276;1218,1221,1223,1224,1225,1227,1228,1231,1233,1235,1236,1238,1243,1248,1253,1257,1262,1267,1272,1277;1222,1225,1226,1227,1229,1230,1231,1234,1236,1238,1239,1241,1246,1251,1255,1260,1264,1269,1274,1279;1226,1228,1230,1231,1232,1233,1235,1237,1239,1241,1242,1244,1249,1253,1258,1262,1266,1271,1276,1280;1229,1232,1233,1234,1235,1237,1238,1240,1242,1244,1245,1247,1251,1256,1260,1264,1268,1273,1278,1282;1233,1235,1236,1238,1239,1240,1241,1243,1246,1247,1248,1250,1254,1258,1262,1266,1270,1275,1279,1284;1236,1239,1240,1241,1242,1243,1244,1246,1249,1250,1251,1253,1257,1261,1265,1269,1272,1277,1281,1285;1240,1242,1243,1244,1245,1246,1247,1249,1252,1253,1254,1256,1260,1264,1267,1271,1275,1279,1283,1287;1243,1245,1246,1247,1248,1249,1250,1252,1255,1256,1256,1258,1262,1266,1270,1273,1277,1281,1285,1289];
    y = y';
    val = interpinbounds(T,P,y,temp,press);
    cp = val;
end

function val = interpinbounds(T,P,y,temp,press)
    % provides 2D interpolation of temperature and pressure
    % with values falling outside of the T and P data being taken as the
    % value at the closest boundary (always returns a valid value)
    val = interp2(T,P,y,temp,press,'linear');
    
    bind = 1*(temp<T(1));
    bind = bind + 2*(temp>T(end));
    bind = bind + 3*(press<P(1));
    bind = bind + 6*(press>P(end));
    
    val(bind==4) = y(end,1);
    val(bind==5) = y(end,end);
    val(bind==7) = y(1,1);
    val(bind==8) = y(1,end);
    val(bind==1) = interp1(P(:,1),y(:,1),press(bind==1));
    val(bind==2) = interp1(P(:,end),y(:,end),press(bind==2));
    val(bind==3) = interp1(T(end,:),y(end,:),temp(bind==3));
    val(bind==6) = interp1(T(1,:),y(1,:),temp(bind==6));
end