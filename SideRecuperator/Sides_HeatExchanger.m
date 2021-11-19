function [Q,m_dot_C,m_dot_H,T_C_out,T_H_out,T_HX,T_Cold,T_Hot,model,results] = Sides_HeatExchanger(opr_dir,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%% Run Specification

% initialize function controls
GeomMode = 'S'; % 'S' or 'H'
BCSpec = 'DP';%'DP' or 'mdot'
fpropmode = 'real';%'isobar' or 'const' or 'real'
spropmode = 'real';%'isobar' or 'const' or 'real'
breakpoints = false; %will pause execution at points, press F5 to continue
plotting = false;
verbose = false;
H_max = 0.01;
maxiters = 5;
fromScratch =  false;
dir_save = opr_dir;
configuration = 'Hthru'; %'Hthru' or 'Cthru'
resistmode = 'Colburn'; %'Colburn';
frictmode = 'Re';%'Re','const'
f_D = 0.181; 
j_C = 8e-3;
global fcorrloc jcorrloc c_kfact_core
fcorrloc = [1,0];
fcorrloc = [0.239087941700866,0.210269930642014];
jcorrloc = [1,0];
replot = false;
results = [];
adaptiveMeshMode = []; %'Heating'

% initialize HX geometric parameters
N_p = 51;  % total number of plates
th_p = 1.5e-3; % thickness of a plate
phi_core = 0.2433; % volume fraction of channels in the core
Wc = 0.2; % width of the core
Lc = 0.6; % length of the core
Ws = 0.006; % thickness of side walls (to which the C side headers are welded)
Ls = 0.006; % thickness of end walls (to which the H side headers are welded)
Hs = 0.006; % thickness of top/bottom plates
D_h = 1.012e-3; % channel hydraulic diameter in m
c_kfact_core = [0.8186, 0.7966, 0.6420]; % global kfact_core

% initialize boundary conditions
m_dot = 0.05;
m_dot_C = m_dot;
m_dot_H = m_dot;
T_C_in = 30; % cold inlet temp in C
T_C_out = 150; % predicted cold outlet temp in C (not needed, will be calculated)
T_H_in = 160; % hot inlet temp in C
T_H_out = 35; % predicted hot outlet temp in C (not needed, will be calculated)
P_C_in = 20e6; % cold inlet pressure in Pa
DP_C = 20e3; % cold pressure drop in Pa
P_C_out = P_C_in-DP_C; % cold outlet pressure in Pa
P_H_in = 8e6; % cold inlet pressure in Pa
DP_H = 20e3; % cold pressure drop in Pa
P_H_out = P_H_in-DP_H; % cold outlet pressure in Pa

%% Parse Function Input

if length(varargin)>0
for i = 1:2:length(varargin)
   switch varargin{i}
       case 'Geometry'
           switch varargin{i+1}
                case {'S','sides'}
                    GeomMode = 'S';
                case {'H','half-sides'}
                    GeomMode = 'H';
            end
            
        case {'Wc','Lc','Ws','Ls','Hs','D_h','N_p','th_p'}
            feval(@()assignin('caller',varargin{i},varargin{i+1}))
        case {'T_C_in','T_C_out','T_H_in','T_H_out','P_C_in','P_C_out','P_H_in','P_H_out','m_dot_C','m_dot_H'}
            feval(@()assignin('caller',varargin{i},varargin{i+1})) 
        case {'DP_C','DP_H'}
%            feval(@()assignin('caller',[varargin{i}(2:4),'_out'],([varargin{i}(2:4),'_in']-varargin{i+1})))
            feval(@()assignin('caller',varargin{i},varargin{i+1}))
        case 'm_dot'
            m_dot = varargin{i+1};
            m_dot_C = m_dot;
            m_dot_H = m_dot;
        case {'H_max'}
           feval(@()assignin('caller',varargin{i},varargin{i+1}))
        case 'MaxIter'
           maxiters = varargin{i+1};
        case {'Replot','replot'}
            switch varargin{i+1}
                case {true,'on','true',1}
                    replot = true;
            end
       case 'flowBC'
           BCSpec = varargin{i+1};
       case 'FluidProperties'
           fpropmode = varargin{i+1};
       case 'SolidProperties'
           spropmode = varargin{i+1};
       case 'breakpoints'
           switch varargin{i+1}
               case {true,'true',1,'on'}
                breakpoints = true;
           end
       case 'plotting'
           switch varargin{i+1}
               case {true,'true',1,'on'}
                plotting = true;
               otherwise
                   plotting = false;
           end
       case 'verbose'
           switch varargin{i+1}
               case {true,'true',1,'on'}
                verbose = true;
           end

       case 'InitialConditions'
           switch varargin{i+1}
               case {'scratch','off'}
                fromScratch = true;
               case {'previous','load','on'}
                fromScratch = true; 
           end
       case 'Nusselt'
           switch varargin{i+1}
               case {true,'true',1,'on'}
                resistmode = 'Nusselt';
           end
       case 'fcorr'
            fcorrloc = varargin{i+1};
       case 'jcorr'
            jcorrloc = varargin{i+1};
       case 'AddaptiveMeshRefinement'
           adaptiveMeshMode = varargin{i+1};
   end
       
end
end

switch BCSpec
    case 'mdot'
    otherwise
        P_H_out = P_H_in-DP_H;
        P_C_out = P_C_in-DP_C;
end

switch GeomMode
case 'H'
m_dot_C = m_dot_C*0.5;
m_dot_H = m_dot_H*0.5;
end


plotting = plotting | replot;

basename = 'HHX_2D_Results';
filename = [dir_save,'\',basename,'_solu.mat'];



warning off

%% Problem Init

N = 5; % default is just the thermal conduction PDE, PDE system size of N=5
model = createHHXT(N);
model.WorkingDirectory = opr_dir;
model.AnalysisType = 'steadystate';
order = 'linear';

%% Geometry
% specify geometry

% determine geometric parameters
N_C = ceil(N_p/2); % number of cold plates
N_H = floor(N_p/2); % number of hot plates
Hc = N_p*th_p; % core height
H = Hc+2*Hs; % total height

% set volume fractions
global phi_C phi_H phi_F phi_S
phi_C = phi_core*(Hc/H)*(N_C/N_p);
phi_H = phi_core*(Hc/H)*(N_H/N_p);
phi_F = phi_C+phi_H;
phi_S = 1-phi_F;

%build the geometry
%and calculate fluid regions, and BC locations
switch GeomMode
case 'H'
[rid_Hh,rid_Ch,rid_Cv,rid_solid,C_in,C_out,H_in,H_out] = HXGeom_halfH(model,Wc,Lc,Ws,Ls);
case 'S'
[rid_Hh,rid_Ch,rid_Cv,rid_solid,C_in,C_out,H_in,H_out] = HXGeom_S(model,Wc,Lc,Ws,Ls);
end

if ~replot
model.setThickness(H);  % set the 2D thickness
if plotting == true
pdegplot(model,'CellLabels','on','EdgeLabels','on','FaceLabels','on','FaceAlpha',0.5);
end
if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end
end

% which stream is hot/cold
switch configuration
    case 'Hthru'
        str_C = 1;
        str_H = 2;
    case 'Cthru'
        str_C = 2;
        str_H = 1;
end
    

%% Solid Material Definition
% define core solid materials

switch spropmode
    case 'const'
        c_k = 17.6; % 316H at 250 C

        c_kfact = [1;0;0;1]; % not using z dimension
        c_k_S = c_kfact*c_k;

        c_kfact = phi_F*c_kfact_core([1,2]) + phi_S*[1,1]; % not using z dimension
        c_kfact = diag(c_kfact); c_kfact = c_kfact(:);
        c_k_ChHh = c_kfact*c_k;

        c_kfact = phi_C*c_kfact_core([2,1]) + phi_H*c_kfact_core([1,2]) + phi_S*[1,1]; % not using z dimension
        c_kfact = diag(c_kfact); c_kfact = c_kfact(:);
        c_k_CvHh = c_kfact*c_k;

        c_kfact = phi_H*c_kfact_core([1,2]) + (phi_S+phi_C)*[1,1]; % not using z dimension
        c_kfact = diag(c_kfact); c_kfact = c_kfact(:);
        c_k_Hh = c_kfact*c_k;
                
        c_kfact = phi_C*c_kfact_core([2,1]) + (phi_S+phi_H)*[1,1]; % not using z dimension
        c_kfact = diag(c_kfact); c_kfact = c_kfact(:);
        c_k_Cv = c_kfact*c_k;
        
    case {'isobar','real'}
        c_k_S = @k_S;
        c_k_ChHh = @k_ChHh;
        c_k_CvHh = @k_CvHh;
        c_k_Hh = @k_Hh;
        c_k_Cv = @k_Cv; 
end
        
% 100% solid regions (no channels)
sld1 = HHXTProperties(model,'Face',rid_solid(1),'Stream',0,...
    'PhaseFraction',1,...
    'ThermalConductivity',c_k_S,...
    'MassDensity',8030,...
    'SpecificHeat',502);
for rid = rid_solid(2:end)
HHXTProperties(model,'Face',rid,'Material',sld1);
end

% regions where both channels exist
rid_counter = intersect(rid_Ch,rid_Hh);
for rid = rid_counter
HHXTProperties(model,'Face',rid,'Material',sld1,...
    'PhaseFraction',phi_S,...
    'ThermalConductivity',c_k_ChHh);
end
rid_cross = intersect(rid_Cv,rid_Hh);
for rid = rid_cross
HHXTProperties(model,'Face',rid,'Material',sld1,...
    'PhaseFraction',phi_S,...
    'ThermalConductivity',c_k_CvHh);
end

% regions where only one channel exists
rid_2chan = intersect([rid_Ch,rid_Cv],rid_Hh);
for rid = rid_Hh(~ismember(rid_Hh,rid_2chan))
HHXTProperties(model,'Face',rid,'Material',sld1,...
    'PhaseFraction',(1-phi_H),...
    'ThermalConductivity',c_k_Hh);
end
for rid = rid_Cv(~ismember(rid_Cv,rid_2chan))
HHXTProperties(model,'Face',rid,'Material',sld1,...
    'PhaseFraction',(1-phi_C),...
    'ThermalConductivity',c_k_Cv);
end

% determine conductive resistance through solid
if isa(c_k_S,'function_handle')
    kmax = c_k_S([],[],20); %max solid conductivity is at 20 C
else
    kmax = c_k_S;
end
kmax = kmax(1);
R_V_cond = 2.8002e-6/kmax; %conductive resistance of the solid
R_V_cond = 0.5*R_V_cond;  %half of resistance goes to each side


%% Fluid Stream Properties
% define fluid stream materials

% use a directional friction factor to keep flow in the horizontal and
% vertical
switch frictmode
    case 'const'
        fwall = f_D*1000000;
        fmatrh = [f_D;fwall;fwall];
        fmatrv = [fwall;f_D;fwall];
    case 'Re'
        fmatrh = @fcorh;
        fmatrv = @fcorv;
end

% how to settup heat transfer solver
switch resistmode
    case 'Colburn'
       rtag = 'jColburn'; 
    case 'Nusselt'
       rtag = 'Nusselt'; 
end

switch fpropmode
    case 'const'
        c_k = 0.055;
        c_rho = 400;
        c_cp = 1600;
        c_mu = 0.000035;
    case 'isobar'
        c_k = @CO2_k;
        c_rho = @CO2_rho;
        c_cp = @CO2_cp;
        c_mu = @CO2_mu;
    case 'real'
        c_k = @CO2_k_real;
        c_rho = @CO2_rho_real;
        c_cp = @CO2_cp_real;
        c_mu = @CO2_mu_real;
end


% in counterflow core

% define cold stream (CO2) material
fldC = HHXTProperties(model,'Face',rid_Ch(1),'Stream',str_C,...
    'PhaseFraction',phi_C,...
    'ThermalConductivity',c_k,...       %thermal conductivity in W/m-K
    'MassDensity',c_rho,...     %densitiy in kg/m3
    'SpecificHeat',c_cp,...      %heat capacity in J/kg-K
    'HydraulicDiameter',D_h,...    %hydraulic diameter in m
    'Viscosity',c_mu,...      %viscosity in kg/m-s
    rtag,@htcoeff_C,...
    'fDarcy',fmatrh,...
    'Resistance',R_V_cond); 
for rid = rid_Ch(2:end)
HHXTProperties(model,'Face',rid,'Material',fldC);
end
fldC = HHXTProperties(model,'Face',rid_Cv(1),'Material',fldC,...
    'fDarcy',fmatrv);
for rid = rid_Cv(2:end)
HHXTProperties(model,'Face',rid,'Material',fldC);
end

% define hot stream (CO2) material
fldH = HHXTProperties(model,'Face',rid_Hh(1),'Stream',str_H,...
    'PhaseFraction',phi_H,...
    'ThermalConductivity',c_k,...       %thermal conductivity in W/m-K
    'MassDensity',c_rho,...     %densitiy in kg/m3
    'SpecificHeat',c_cp,...      %heat capacity in J/kg-K
    'HydraulicDiameter',D_h,...    %hydraulic diameter in m
    'Viscosity',c_mu,...      %viscosity in kg/m-s
    rtag,@htcoeff_H,...
    'fDarcy',fmatrh,...
    'Resistance',R_V_cond);
for rid = rid_Hh(2:end)
HHXTProperties(model,'Face',rid,'Material',fldH);
end

%% Mesh and Prepare model

% generate a mesh of tets
if verbose == true
disp('    ...generating mesh')
end
% generateMesh(model,'GeometricOrder',order);
generateMesh(model,'Hmax',H_max,'GeometricOrder',order);

if true
for i = 1
    elem_refine = [];
    for EdgeID = [C_in,C_out,H_in,H_out]
        nodes = model.Mesh.findNodes('region','Edge',EdgeID);
        elem_edge = model.Mesh.findElements('attached',nodes);
        elem_refine = [elem_refine;elem_edge']; % need to be passed as a column
    end
    model = refineHHXTmesh(model,elem_refine);  % elements to refine, input as a column
end
end

if plotting == true
pdeplot(model,'FaceAlpha',0.5);
end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

%% Apply Initial Conditions
% set initial conditions
% these are just guesses so the solver can start from somewhere
% we will guess average values
u0 = [0.25*(T_C_in+T_C_out+T_H_in+T_H_out),...
      0.5*(P_C_in+P_C_out),0.5*(T_C_in+T_C_out),...
      0.5*(P_H_in+P_H_out),0.5*(T_H_in+T_H_out)]';
setInitialConditions(model,u0);

%% Load Previous Solution or Perform Init Step

% see if previous solution exists and if it has the same mesh
test = true;
if isfile([filename,'.mat'])
    load(filename,'results');
    test = size(results.Mesh.Nodes,2) ~= size(model.Mesh.Nodes,2);    
end

if fromScratch || test % run an init step: simple conductive PDE with Dirichlet BCs

    %boundary conditions     
    T_sep = 5;
    DBC = [T_C_in,P_C_in,T_C_in-T_sep];
    EBC = [1,2,3];
    applyBoundaryCondition(model,'mixed','Edge',C_in,'u',DBC,'EquationIndex',EBC);
    DBC = [T_H_out,P_H_out,T_H_out+T_sep];
    EBC = [1,4,5];
    applyBoundaryCondition(model,'mixed','Edge',H_out,'u',DBC,'EquationIndex',EBC);
    DBC = [T_H_in,P_H_in,T_H_in+T_sep];
    EBC = [1,4,5];
    applyBoundaryCondition(model,'mixed','Edge',H_in,'u',DBC,'EquationIndex',EBC);
    DBC = [T_C_out,P_C_out,T_C_out-T_sep];
    EBC = [1,2,3];
    applyBoundaryCondition(model,'mixed','Edge',C_out,'u',DBC,'EquationIndex',EBC);

    if verbose == true
    disp('starting initialization: from simple conductive model...')  
    end

    specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',zeros(N,1));

    % solve step    
    if verbose == true
    tic;
    disp('    ...solving step')
    end
    results = solvepde(model);
    if verbose == true
    solve_time = toc;
    finish_date = datestr(now);
    disp(['finished initialization: at ',finish_date,' elapsed time ',num2str(solve_time),' sec'])
    end

else % Load previously saved solution    
    if verbose == true
    disp('starting initialization: loading from solution...')  
    end
    load(filename,'results')
    if verbose == true
    disp(['finished initialization: loaded ',filename,'.mat'])
    end

end

%% Define Stream Boundary Conditions
delete(model.BoundaryConditions);

%boundary conditions       
switch BCSpec   
    case 'DP'
        %cold stream
        applyStreamBoundaryCondition(model,'Edge',C_in,'Stream',str_C,'Direction','inlet',...
            'Pressure',P_C_in,'Temperature',T_C_in);
        applyStreamBoundaryCondition(model,'Edge',C_out,'Stream',str_C,'Direction','outlet',...
            'Pressure',P_C_out);
        %hot stream
        applyStreamBoundaryCondition(model,'Edge',H_in,'Stream',str_H,'Direction','inlet',...
            'Pressure',P_H_in,'Temperature',T_H_in);
        applyStreamBoundaryCondition(model,'Edge',H_out,'Stream',str_H,'Direction','outlet',...
            'Pressure',P_H_out);    
    case 'mdot'
        %cold stream
        applyStreamBoundaryCondition(model,'Edge',C_in,'Stream',str_C,'Direction','inlet',...
            'Pressure',P_C_in,'Temperature',T_C_in,'MassFlow',m_dot_C);
        applyStreamBoundaryCondition(model,'Edge',C_out,'Stream',str_C,'Direction','outlet',...
            'MassFlow',m_dot_C);
        %hot stream
        applyStreamBoundaryCondition(model,'Edge',H_in,'Stream',str_H,'Direction','inlet',...
            'Pressure',P_H_in,'Temperature',T_H_in,'MassFlow',m_dot_H);
        applyStreamBoundaryCondition(model,'Edge',H_out,'Stream',str_H,'Direction','outlet',...
            'MassFlow',m_dot_H);       
end


%% Steady State Solution: coupled streams with variable properties and calculated permeability

%% Steady State Step 1: coupled streams with variable properties and calculated permeability
if verbose == true
disp('starting steady state solution :...')    
end

if ~isempty(results)
setInitialConditions(model,results);
end

% solver options
options = {'Permeability','calculated','ResistanceMode','mixed', ...
           'Conductivity','aniso','FluidConduction','on','Advection','on'};
if ~isempty(adaptiveMeshMode) 
options{length(options)+1} = 'AddaptiveMeshRefinement';
% add max of 2 steps of Addaptive mesh refinement for Heating (q_dot)
options{length(options)+1} = {2,'Heating'}; 
end
model.MatrixOptions = options;

if verbose == true             
    model.SolverOptions.ReportStatistics = 'on';
else
   model.SolverOptions.ReportStatistics = 'off'; 
end
model.SolverOptions.MaxIterations = maxiters;
model.SolverOptions.ResidualTolerance = 1e-6;
model.SolverOptions.MinStep = 0.5; % smaller steps (<1) can be taken but are not necessary
model.SolverOptions.AbsoluteTolerance = 1e-2;
model.SolverOptions.ResidualNorm = 'energy';

if verbose == true
tic;
disp('    ...solving step')
end
results = solveHHXTpde(model);
if verbose == true
solve_time = toc;
finish_date = datestr(now);
disp(['finished steady state solution : at ',finish_date,' elapsed time ',num2str(solve_time),' sec'])
end

% determine velocities and mass fluxes within the core
vel_interps = PullVelocities(model,results,0.5*Lc,0.5*Wc);
mdot_interps = PullMassFlows(model,results,c_rho,0.5*Lc,0.5*Wc);

% save(filename,'model','results','vel_interps','mdot_interps');
save(filename,'model','results','vel_interps','mdot_interps','-nocompression');

if verbose == true
disp(['    state saved to ',filename,'.mat']);
end

%% Display Results

switch configuration
    case 'Hthru'
        str_C = 1;
        str_H = 2;
    case 'Cthru'
        str_C = 2;
        str_H = 1;
end

if verbose == true
% display results
disp('Results summary:')
disp(results.Tables{1})

disp('')
disp('Cold Stream:')
disp(results.Tables{str_C+2})

disp('')
disp('Hot Stream:')
disp(results.Tables{str_H+2})

disp('')
disp('Number of Elements:')
disp(length(results.Mesh.Elements))
end

if plotting == true
PlotResults(filename,dir_save,breakpoints,R_V_cond,[str_C,str_H]);
PlotMassFlowDists(filename,dir_save,breakpoints,c_rho,[str_C,str_H],rid_counter);
PlotVelocityDists(filename,dir_save,breakpoints,[str_C,str_H],rid_counter);
end

warning on

%% Output

T_HX = results.Tables{1};
T_Cold = results.Tables{str_C+2};
T_Hot = results.Tables{str_H+2};

switch GeomMode
    case 'H'
Q = 2*mean(abs(T_HX.Heating)); % total heat transfered in W (model is a half symmetry so x2)
m_dot_C = 2*mean(abs(T_Cold.Mass_Flow)); % cold side mass flow in kg/s (model is a half symmetry so x2)
m_dot_H = 2*mean(abs(T_Hot.Mass_Flow)); % hot side mass flow in kg/s (model is a half symmetry so x2)
T_C_out = T_Cold.Temperature(2); % cold outlet temp as calculated by model
T_H_out = T_Hot.Temperature(2); % hot outlet temp as calculated by model
    case 'S'
Q = mean(abs(T_HX.Heating)); % total heat transfered in W
m_dot_C = mean(abs(T_Cold.Mass_Flow)); % cold side mass flow in kg/s 
m_dot_H = mean(abs(T_Hot.Mass_Flow)); % hot side mass flow in kg/s
T_C_out = T_Cold.Temperature(2); % cold outlet temp as calculated by model
T_H_out = T_Hot.Temperature(2); % hot outlet temp as calculated by model        
end

%% Conductivity Functions

function k = k_matl(temp)
    k = M316_k_real([],[],temp);
end

function kmatrix = k_S(pos,time,temp)
    kfact = [1;0;0;1]; % not using z dimension
    k = k_matl(temp);
    kmatrix = kfact*k;
end

function kmatrix = k_ChHh(pos,time,temp)
%     global phi_F phi_S
    kfact = phi_F*c_kfact_core([1,2]) + phi_S*[1,1]; % not using z dimension
    kfact = diag(kfact); kfact = kfact(:);
    k = k_matl(temp);
    kmatrix = kfact*k;
end

function kmatrix = k_CvHh(pos,time,temp)
%     global phi_C phi_H phi_S
    kfact = phi_C*c_kfact_core([2,1]) + phi_H*c_kfact_core([1,2]) + phi_S*[1,1]; % not using z dimension
    kfact = diag(kfact); kfact = kfact(:);
    k = k_matl(temp);
    kmatrix = kfact*k;
end

function kmatrix = k_Hh(pos,time,temp)
%     global phi_C phi_H phi_S
    kfact = phi_H*c_kfact_core([1,2]) + (phi_S+phi_C)*[1,1]; % not using z dimension
    kfact = diag(kfact); kfact = kfact(:);
    k = k_matl(temp);
    kmatrix = kfact*k;
end

function kmatrix = k_Cv(pos,time,temp)
%     global phi_C phi_H phi_S
    kfact = phi_C*c_kfact_core([2,1]) + (phi_S+phi_H)*[1,1]; % not using z dimension
    kfact = diag(kfact); kfact = kfact(:);
    k = k_matl(temp);
    kmatrix = kfact*k;
end

%% Friction Factor Functions

function fmatr = fcorh(Re)
    f = f_ANL(Re(1,:));
    f = f*fcorrloc(1).*Re(1,:).^fcorrloc(2);
    fmatr = [f;f*1000000];
end

function fmatr = fcorv(Re)
    f = f_ANL(Re(2,:));
    f = f*fcorrloc(1).*Re(2,:).^fcorrloc(2);
    fmatr = [f*1000000;f];
end

%% Heat Transferm Functions

function htcoeff = htcoeff_C(pos,time,Re,Pr)
    htcoeff = j_UWC(Re,Pr);
    htcoeff = htcoeff*jcorrloc(1).*Re.^jcorrloc(2);
end

function htcoeff = htcoeff_H(pos,time,Re,Pr)
    htcoeff = j_UWH(Re,Pr);
    htcoeff = htcoeff*jcorrloc(1).*Re.^jcorrloc(2);
end

end

