function [f_D_c,j_C_c] = FJ_optimizationSearch(dir_save,HX_state,search_mode,f_D_lims,j_C_lims,f_D_start,j_C_start,useNu)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%%!!!! start with low values on these !!!
options = optimset('MaxIter',20);
max_iter = 1;
% options = optimset('MaxIter',500);
% max_iter = 10;

if isempty(search_mode)
search_mode = 'fiber-hydro'; % 'temp-hydro-qdot' %options are: fiber, hydro, temp, qdot (seperate by -)
end
split = strsplit(search_mode,'-');

% global f_D_c j_C_c
f_D_c = f_D_start;
j_C_c = j_C_start;  %C_Nu ~ 17x j_C

% save('fjCalc','f_D_c','j_C_c');

% options = optimset('Display','iter');
% options = optimset('UseParallel',true);
% options = optimoptions('solvername','UseParallel',true);

% options = optimset('UseParallel',true,'MaxIter',10);

% options = optimoptions('fmincon','UseParallel',true,'MaxIter',10);

[~,results] = ZZ_Recuperator(dir_save,'HX_state',HX_state,...\
    'j_C',j_C_c,'f_D',f_D_c,'CO2','real','MaxIter',15,'Nusselt',useNu,...
    'InitialConditions','scratch');

% disp('starting optimization: ...')  
% tic;
f_D_last = f_D_c;
j_C_last = j_C_c;
tol = 1e-3;

conv = false;
iter = 1;
while ~conv
    
    for i = 1:length(split)
        switch split{i}
            case 'fiber'
                % search j_D to minimize outlet fiber temperature deviation
%                 j_C_c = fminbnd(@minimize_fiberdiff,j_C_lims(1),j_C_lims(2),options);
                j_C_c = fminsearch(@minimize_fiberdiff,j_C_c,options);
                % save('fjCalc','f_D_c','j_C_c');
                
            case 'temp'
                % search j_D to minimize outlet temperature deviation
%                 j_C_c = fminbnd(@minimize_tempdiff,j_C_lims(1),j_C_lims(2),options);
                j_C_c = fminsearch(@minimize_tempdiff,j_C_c,options);
                % save('fjCalc','f_D_c','j_C_c');
                
            case 'hydro'
                % search f_C to minimize pressure drop deviation
%                 f_D_c = fminbnd(@minimize_DPdiff,f_D_lims(1),f_D_lims(2),options);
                f_D_c = fminsearch(@minimize_DPdiff,f_D_c,options);
                % save('fjCalc','f_D_c','j_C_c');
                
            case 'qdot'
                % search j_D to minimize q_dot deviation
%                 j_C_c = fminbnd(@minimize_qdotdiff,j_C_lims(1),j_C_lims(2),options);
                j_C_c = fminsearch(@minimize_qdotdiff,j_C_c,options);
                % save('fjCalc','f_D_c','j_C_c');
                
        end
    end
    
    j_C_c = abs(j_C_c);
    f_D_c = abs(f_D_c);

    % test to see if j and f have converged
    test = [abs((f_D_c-f_D_last)/f_D_last),abs((j_C_c-j_C_last)/j_C_last)];
    conv = (max(test) <= tol) || (iter >= max_iter);

    f_D_last = f_D_c;
    j_C_last = j_C_c;
    iter = iter +1;
end

%%report results

% solve_time = toc;
% finish_date = datestr(now);
% disp(['finished : at ',finish_date,' elapsed time ',num2str(solve_time),' sec'])


% [cmp_fiber(i,:),cmp_TC(i,:),fpnt(i,:),epnt(i,:)] = ModelVSFibers(HX_state,results);
% [cmp_hyd(i,:),cmp_tmp(i,:),cmp_hxt(i,:),cmp_rey(i,:),...
%     hyd(i,:),tmp(i,:),hxt(i,:),rey(i,:)] = ModelVSConditions(HX_state,results,'mdot');

function mintar = minimize_fiberdiff(j_C)
%     global f_D_c
%     load('fjCalc','f_D_c');
    [~,results] = ZZ_Recuperator(dir_save,'HX_state',HX_state,'j_C',abs(j_C),'f_D',f_D_c,'CO2','real','MaxIter',15,'Nusselt',useNu);

    [cmp_fiber,~,~,~] = ModelVSFibers(HX_state,results);
    
    mintar = cmp_fiber(1);

end

function mintar = minimize_tempdiff(j_C)
%     global f_D_c
%     load('fjCalc','f_D_c');
    
    [~,results] = ZZ_Recuperator(dir_save,'HX_state',HX_state,'j_C',abs(j_C),'f_D',f_D_c,'CO2','real','MaxIter',15,'Nusselt',useNu);

    [~,cmp_tmp,~,~,...
        ~,~,~,~] = ModelVSConditions(HX_state,results,'mdot');
    
    mintar = cmp_tmp(1);

end

function mintar = minimize_DPdiff(f_D)
%     global j_C_c
%     load('fjCalc','j_C_c');
    [~,results] = ZZ_Recuperator(dir_save,'HX_state',HX_state,'j_C',j_C_c,'f_D',abs(f_D),'CO2','real','MaxIter',15,'Nusselt',useNu);

    [cmp_hyd,~,~,~,...
        ~,~,~,~] = ModelVSConditions(HX_state,results,'mdot');
    
    mintar = cmp_hyd(1);

end

function mintar = minimize_qdotdiff(j_C)
%     global f_D_c
%     load('fjCalc','f_D_c');
    [~,results] = ZZ_Recuperator(dir_save,'HX_state',HX_state,'j_C',abs(j_C),'f_D',f_D_c,'CO2','real','MaxIter',15,'Nusselt',useNu);

    [~,~,cmp_hxt,~,...
        ~,~,~,~] = ModelVSConditions(HX_state,results,'mdot');
    
    mintar = cmp_hxt(1);

end

end

