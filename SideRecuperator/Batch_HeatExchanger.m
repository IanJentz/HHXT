
% need to have functions folder loaded into the Matlab path
% this allows CO2 propery functions, heat transfer and friction factor
% coefficients, 316 stainless property functions, to be loaded
addpath(genpath([cd,'\functions']))


% choose where to save this run
dir_save = [cd,'\Results'];   

% how model runs is setup using input arguments which are passed in
b_modelargin = {'Geometry','S'}; % choose which geometry to use 'S' is one side entrance S-shape, 'H' is two sides entrance H-shape

% you can keep building up the input argument cell array
b_modelargin = [b_modelargin,{'H_max',0.01,'MaxIter',5}]; % Mesh size (lower -> more elements), and maximum iterations (convergence can happen in less iterations)
b_modelargin = [b_modelargin,{'plotting','on','verbose','off'}]; % turn plotting and display of progress on
    
% set boundary conditions
b_modelargin = [b_modelargin,{'T_C_in',40,'T_H_in',200}]; % inlet temperatures in C
b_modelargin = [b_modelargin,{'T_C_out',190,'T_H_out',45}]; % guess out temperatures in C, for initial conditions
b_modelargin = [b_modelargin,{'P_C_in',20e6,'P_H_in',10e6}]; % inlet pressures in Pa
m_dot = 0.4; % mass flow in kg/s


L_var = 0.2:0.1:0.8;  % vary the core length (in m) of the HX
W_var = 0.05:0.05:0.2;  % vary the core width (in m) of the HX
Np_var = [31,51,71];  % vary the number of plates stacked in the HX (only use odd numbers)

[L_var,W_var,Np_var]  = meshgrid(L_var,W_var,Np_var);
L_var = L_var(:); W_var = W_var(:); Np_var = Np_var(:);
rmind = 2*W_var > L_var; % if W > 1/2, the side entrances won't fit, these are invalide geometries
L_var(rmind) = []; W_var(rmind) = []; Np_var(rmind) = [];
N_tr = length(L_var); % number of trials to run

% pre-allocate a Table that will be filled out in the par-for loop
varnames = {'Q','m_dot_C','m_dot_H','T_C_out','T_H_out'};
vartypes = varnames; for i = 1:length(vartypes); vartypes{i} = 'double'; end
T_results = table('Size',[N_tr,length(varnames)],'VariableTypes',vartypes,'VariableNames',varnames);

% parfor loops need arrays to be passed in and out of the loop
% these have to by intially sized
ME_c = cell(N_tr,1);
p_Q = zeros(N_tr,1);
p_m_dot_C = p_Q; p_m_dot_H = p_Q; p_T_C_out = p_Q; p_T_H_out = p_Q;

tic;
parfor i = 1:N_tr
dir_opr = [dir_save,'\trial_',num2str(i,'%03.0f')];
if ~isdir(dir_opr); mkdir(dir_opr); end
    
% heat exchanger size
i_modelargin = [b_modelargin,{'Lc',L_var(i),'Wc',W_var(i)}]; % length and width of core
i_modelargin = [i_modelargin,{'N_p',Np_var(i),'th_p',1.5e-3}]; % number of plates stacked (C-H-C-H ... C-H-C) and plate thickness -> results in total height
i_modelargin = [i_modelargin,{'Ls',0.01,'Ws',0.01,'Hs',0.01}]; % side wall, end wall, and top/bottom plate thickness

% first run with constant mass flow boundary conditions
i_modelargin = [i_modelargin,{'InitialConditions','scratch'}]; % start solution from scratch
i_modelargin = [i_modelargin,{'flowBC','mdot'}]; % start with constant mass flux boundary conditions
i_modelargin = [i_modelargin,{'m_dot_C',m_dot,'m_dot_H',m_dot}]; % set mass flow values

try

[~,~,~,~,~,~,T_Cold,T_Hot,~,~] = ...
    Sides_HeatExchanger(dir_opr,i_modelargin{:}); % build and solve model
% determine the pressure drops
DP_C = abs(diff(T_Cold.Pressure));
DP_H = abs(diff(T_Hot.Pressure));


% finally run with uniform pressure drop, the natural boundary condition
i_modelargin = [i_modelargin,{'InitialConditions','previous'}]; % use previous solution as IC
i_modelargin = [i_modelargin,{'flowBC','DP'}]; % swithc to natural DP boundary condition
i_modelargin = [i_modelargin,{'DP_C',DP_C,'DP_H',DP_H}]; % pressure drop in Pa

[p_Q(i),p_m_dot_C(i),p_m_dot_H(i),p_T_C_out(i),p_T_H_out(i),T_HX,T_Cold,T_Hot,model,results] = ...
    Sides_HeatExchanger(dir_opr,i_modelargin{:}); % build and solve model

% update progress
disp(['finished trial ',num2str(i),' of ',num2str(N_tr)])

catch ME
    ME_c{i} = ME;
    disp(['failed at trial ',num2str(i),' of ',num2str(N_tr)])
    failed_runs(i) = true;
end

end
toc;
solve_time = toc;
finish_date = datestr(now);
disp(['finished trials : at ',finish_date,' elapsed time ',num2str(solve_time),' sec'])

% add results to table
T_results.Q = p_Q;
T_results.m_dot_C = p_m_dot_C;
T_results.m_dot_H = p_m_dot_H;
T_results.T_C_out = p_T_C_out;
T_results.T_H_out = p_T_H_out;





