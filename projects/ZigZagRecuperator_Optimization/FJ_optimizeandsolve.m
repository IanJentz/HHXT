function [model,results,Re,fj] = FJ_optimizeandsolve(HX_state,recupe_type,search_mode,dir_save,useNu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    verbose = false;
    dir_save =cd;
else
    verbose = false;
    warning off
    mkdir([dir_save,'/figs']);
    warning on
end


S = load('C:\Users\ANSYS_Solver\Documents\AFandZZ\postProcessing\OptimizationSpace.mat');

% switch useNu
%     case 'on'
%         S.guess(2) = 15*S.guess(2);
%         S.bounds(3:4) = (10^0.5)*1e1*S.bounds(3:4); 
% end

[f_D_c,j_C_c] = FJ_optimizationSearch(dir_save,HX_state,search_mode,S.bounds(1:2),S.bounds(3:4),S.guess(1),S.guess(2),useNu);

switch recupe_type
    case 'ZZ'
[model,results] = ZZ_Recuperator(dir_save,...
    'HX_state',HX_state,'f_D',f_D_c,'j_C',j_C_c,...
    'CO2','real','MaxIter',100,...
    'verbose',verbose,'plotting','on',...
    'H_max',0.005,'Nusselt',useNu);
    case 'AF'
        error('airfoil model not yet defined');
end
[cmp_fiber,cmp_TC,~,~] = ModelVSFibers(HX_state,results,dir_save);

Re = [results.Tables{1}.Reynolds(1),results.Tables{1}.Reynolds(2)];
fj = [f_D_c,j_C_c];
S.guess = fj;
% save('W:\PCHE\AFandZZ\postProcessing\OptimizationSpace.mat','-struct','S')


end

