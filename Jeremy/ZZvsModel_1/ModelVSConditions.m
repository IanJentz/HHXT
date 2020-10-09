function [cmp_hyd,cmp_tmp,cmp_hxt,cmp_rey,hyd_model,tmp_model,hxt_model,rey_model] = ModelVSConditions(obj,results,BCmode)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

BC_ind = 3:4;
BCloc = 'outlet';
rout = 0*BC_ind;
for i = 1:length(BC_ind)
    rout(i) = find(strcmp(results.Tables{1,BC_ind(i)}.Direction,'outlet'));
end
    
hyd_model = 0*BC_ind; 
tmp_model = 0*BC_ind; 
switch BCmode
    case 'mdot'
        for i = 1:length(BC_ind)
            hyd_model(i) = results.Tables{1,1}.Pressure_Drop(i);
%             hyd_model(i) = results.Tables{1,BC_ind(i)}.Pressure(rout(i));
            tmp_model(i) = results.Tables{1,BC_ind(i)}.Temperature(rout(i));
        end
        hyd_state = [ obj.DP_C.mean ,obj.DP_H.mean ];
        if strcmp(obj.DP_C.unit,'psi')
            hyd_state = 6894.76*hyd_state;
        end
        tmp_state = [ obj.T_C_out.mean , obj.T_H_out.mean ];
            
    case 'DP'
        for i = 1:length(BC_ind)
            hyd_model(i) = results.Tables{1,BC_ind(i)}.Mass_Flow(rout(i));
            tmp_model(i) = results.Tables{1,BC_ind(i)}.Temperature(rout(i));
        end
        hyd_state = [ obj.m_dot_C.mean , obj.m_dot_H.mean];
        if strcmp(obj.P_C_in.unit,'kg/h')
            hyd_state = 3600*hyd_state;
        end
        tmp_state = [ obj.T_C_out.mean , obj.T_H_out.mean ];
end

hxt_model = 0*BC_ind;
rey_model = 0*BC_ind;
for i = 1:length(BC_ind)
    hxt_model(i) = results.Tables{1,1}.Heating(i);
    rey_model(i) = results.Tables{1,1}.Reynolds(i);
end
if ~isempty(obj.EES_sol)
hxt_state = [ obj.EES_sol.q_dot_C , obj.EES_sol.q_dot_H ];
rey_state = [ obj.EES_sol.Re_bar_C , obj.EES_sol.Re_bar_H ];
else
    hxt_state = [NaN,NaN];
    rey_state = [NaN,NaN];
end

hyd_dif = hyd_model-hyd_state; hyd_dif = abs(hyd_dif);
cmp_hyd = [rms(hyd_dif),std(hyd_dif)];
tmp_dif = tmp_model-tmp_state; tmp_dif = abs(tmp_dif);
cmp_tmp = [rms(tmp_dif),std(tmp_dif)];
hxt_dif = hxt_model-hxt_state; hxt_dif = abs(hxt_dif);
cmp_hxt = [rms(hxt_dif),std(hxt_dif)];
rey_dif = rey_model-rey_state; rey_dif = abs(rey_dif);
cmp_rey = [rms(rey_dif),std(rey_dif)];



end

