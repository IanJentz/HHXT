function [rmsd_fiber,rmsd_TC,rmsd_DP,rmsd_Tout,rmsd_qdot,Re_C,Re_H] = ZZvsModel(j_C,f_D)
%ZZVSMODEL runs the ZZ_Recuperator.m model and compares the result to
%          experimental data given in HX_state.mat.
%   Inputs:
%     j_C - Colburn heat transfer coefficient
%     f_D - Darcy friction factor coefficient
%
%   Outputs:
%     rmsd_fiber - rms difference between the model and fiber readings
%     rmsd_TC    - rms difference between the model and embedded thermocouple readings
%     rmsd_DP    - rms difference between the model and measured pressure drop
%     rmsd_Tout  - rms difference between the model and measured outlet temperatures
%     rmsd_qdot  - rms difference between the model and measured heat transfer
%     Re_C       - cold side Reynolds number
%     Re_H       - hot side Reynolds number


%% setup the problem
% here we choice whether to display progress, plot results, and we load in
% the experimental data HX_state.mat which will be used in the model.

verbose = false;  % choice of whether to print progress to the command line
                  % if would be usefull to have each iteration create a log
                  % which would contain the matlab command line output
plotting = false; % choice of whether to plot results to a figs folder
                  % you might need to get each node to be able to dsiplay
                  % and write graphics in order to use this
if plotting
    warning off
    mkdir([dir_save,'/figs']);
    warning on
end

dir_save = cd; % the save directory might have to be changed for each node 

data = load('HX_state.mat'); % load the experimental data from file
HX_state = data.HX_state; 

%% run the model
% runs the ZZ_Recuperator model using the j_C and f_D that are passed in.
% We use a high refinement mesh with max elem size of 0.005, use real gas
% CO2 properties, and allow the solution to attempt at most 100 iterations.
[model,results] = ZZ_Recuperator(dir_save,...
    'HX_state',HX_state,'f_D',f_D,'j_C',j_C,...
    'CO2','real','MaxIter',100,...
    'verbose',verbose,'plotting',plotting,...
    'H_max',0.005);

%% calculate results
% here we will run some scripts that calculate the model vs. experiment
% metrics which we would like to output

Re_C = results.Tables{1}.Reynolds(1); % cold side Reynolds number
Re_H = results.Tables{1}.Reynolds(2); % hot side Reynolds number

[cmp_fiber,cmp_TC,~,~] = ModelVSFibers(HX_state,results);
rmsd_fiber = cmp_fiber(1); % rms difference between the model and fiber readings
rmsd_TC = cmp_TC(1); % rms difference between the model and embedded thermocouple readings

[cmp_hyd,cmp_tmp,cmp_hxt,~,~,~,~,~] = ModelVSConditions(HX_state,results,'mdot');
rmsd_DP = cmp_hyd(1); % rms difference between the model and measured pressure drop
rmsd_Tout = cmp_tmp(1); % rms difference between the model and measured outlet temperatures
rmsd_qdot = cmp_hxt(1); % rms difference between the model and measured heat transfer

end

