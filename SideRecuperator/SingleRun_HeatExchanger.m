
% need to have functions folder loaded into the Matlab path
% this allows CO2 propery functions, heat transfer and friction factor
% coefficients, 316 stainless property functions, to be loaded
addpath(genpath([cd,'\functions']))


% choose where to save this run
dir_save = [cd,'\Results'];   

% how model runs is setup using input arguments which are passed in
modelargin = {'Geometry','S'}; % choose which geometry to use 'S' is one side entrance S-shape, 'H' is two sides entrance H-shape

% you can keep building up the input argument cell array
modelargin = [modelargin,{'H_max',0.01,'MaxIter',5}]; % Mesh size (lower -> more elements), and maximum iterations (convergence can happen in less iterations)
modelargin = [modelargin,{'plotting','on','verbose','on'}]; % turn plotting and display of progress on

% heat exchanger size
modelargin = [modelargin,{'Lc',0.4,'Wc',0.1}]; % length and width of core
modelargin = [modelargin,{'N_p',31,'th_p',1.5e-3}]; % number of plates stacked (C-H-C-H ... C-H-C) and plate thickness -> results in total height
modelargin = [modelargin,{'Ls',0.01,'Ws',0.01,'Hs',0.01}]; % side wall, end wall, and top/bottom plate thickness
    
% set boundary conditions
modelargin = [modelargin,{'T_C_in',40,'T_H_in',200}]; % inlet temperatures in C
modelargin = [modelargin,{'T_C_out',190,'T_H_out',45}]; % guess out temperatures in C, for initial conditions
modelargin = [modelargin,{'P_C_in',20e6,'P_H_in',10e6}]; % inlet pressures in Pa




% first run with constant mass flow boundary conditions
modelargin = [modelargin,{'InitialConditions','scratch'}]; % start solution from scratch
modelargin = [modelargin,{'flowBC','mdot'}]; % start with constant mass flux boundary conditions
m_dot = 0.4; % mass flow in kg/s
modelargin = [modelargin,{'m_dot_C',m_dot,'m_dot_H',m_dot}]; % set mass flow values

[~,~,~,~,~,~,T_Cold,T_Hot,~,~] = ...
    Sides_HeatExchanger(dir_save,modelargin{:}); % build and solve model
% determine the pressure drops
DP_C = abs(diff(T_Cold.Pressure));
DP_H = abs(diff(T_Hot.Pressure));


% finally run with uniform pressure drop, the natural boundary condition
modelargin = [modelargin,{'InitialConditions','previous'}]; % use previous solution as IC
modelargin = [modelargin,{'flowBC','DP'}]; % swithc to natural DP boundary condition
modelargin = [modelargin,{'DP_C',DP_C,'DP_H',DP_H}]; % pressure drop in Pa

[Q,m_dot_C,m_dot_H,T_C_out,T_H_out,T_HX,T_Cold,T_Hot,model,results] = ...
    Sides_HeatExchanger(dir_save,modelargin{:}); % build and solve model






