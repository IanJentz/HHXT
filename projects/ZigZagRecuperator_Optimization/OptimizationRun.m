
% load and set the starting guesses for f and j
S = load('OptimizationSpace.mat');
S.guess(1) = 0.6; %starting f_D guess
S.guess(2) = 0.02; % starting j guess
save('OptimizationSpace.mat','-struct','S');

% load the experimental state HX_state
load('HX_state.mat')

search_mode = 'fiber-hydro'; %'fiber-hydro'; %options are: fiber, hydro, temp, qdot (seperate by -)
recupe_type = 'ZZ';
dir_save = cd;
useNu = false; % it will run using a constant colburn factor for j, and not a Nusselt number function

[model,results,Re,fj] = FJ_optimizeandsolve(HX_state,recupe_type,search_mode,dir_save,useNu);