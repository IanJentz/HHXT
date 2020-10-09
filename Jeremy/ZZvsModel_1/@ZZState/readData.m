function obj = readData(obj)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here
    addpath('W:\PCHE\AFandZZ\RTDCalibration\TDMSReader')
%     addpath('W:\PCHE\AFandZZ\RTDCalibration\TDMSReader\private')
    addpath('W:\PCHE\AFandZZ\RTDCalibration\TDMSReader\tdmsSubfunctions')

    tdms_struct = TDMS_getStruct(obj.fileData);

    obj.time_start = datetime(tdms_struct.Temperatures_Inst_Times.Time_Stamp.data(1),'ConvertFrom', 'datenum');
    obj.seconds = tdms_struct.Temperatures_Inst_Times.Seconds.data;
    obj.ntimes = length(obj.seconds);

    group = 'Temperatures_Inst';
    obj.T_H_in = TDMSvar(tdms_struct,group,'TC20');    
    obj.T_H_out = TDMSvar(tdms_struct,group,'TC19');
    obj.T_C_in = TDMSvar(tdms_struct,group,'TC22');
    obj.T_C_out = TDMSvar(tdms_struct,group,'TC21');

    WTCs = {'TC48','TC47','TC46','TC45','TC41','TC42','TC43','TC44'};
    obj.pos_W = (25.4/1000)*...
        [ 2.2398 2.2398 2.2398 2.2398 20.335 20.335 20.335 20.335;...
          4.1575 2.9693 1.7811 0.5930 4.1575 2.9693 1.7811 0.5930 ];
    y_T = obj.pos_W(2,1:(size(obj.pos_W,2)/2));
    obj.T_W = cell(length(WTCs),1);

    for i = 1:length(WTCs)
        obj.T_W{i} = TDMSvar(tdms_struct,group,WTCs{i});
    end
    
    obj.pos_TF = (25.4/1000)*[...
        [ 2.2398 2.2398 2.2398 2.2398 2.2398 2.2398 2.2398 2.2398;...
          4.4550 4.1109 3.2668 2.6727 2.0786 1.4845 0.8904 0.2963 ],...
        [ 20.335 20.335 20.335 20.335 20.335 20.335 20.335 20.335;...
          4.4550 4.1109 3.2668 2.6727 2.0786 1.4845 0.8904 0.2963 ] ];
    y_ind1 = 1:(size(obj.pos_W,2)/2);
    y_ind2 = ((size(obj.pos_W,2)/2)+1):size(obj.pos_W,2);
    y_F = obj.pos_TF(2,(1:(size(obj.pos_TF,2)/2)));
    T_F = zeros(size(obj.pos_TF,2),obj.ntimes);
    T_F_uncr = T_F;
    obj.T_F = cell(size(obj.pos_TF,2),1);
    TWm = zeros(length(WTCs),obj.ntimes);
    TWm_uncr = TWm;
    for i = 1:length(WTCs)
        TWm(i,:) = obj.T_W{i}.data;
        TWm_uncr(i,:) = obj.T_W{i}.uncr;
    end
    for j = 1:obj.ntimes
        T_F1 = interp1(y_T,TWm(y_ind1,j)',y_F,'linear','extrap');
        T_F2 = interp1(y_T,TWm(y_ind2,j)',y_F,'linear','extrap');
        T_F(:,j) = [T_F1';T_F2'];
        T_F1_uncr = interp1(y_T,TWm_uncr(y_ind1,j)',y_F,'linear','extrap');
        T_F2_uncr = interp1(y_T,TWm_uncr(y_ind2,j)',y_F,'linear','extrap');
        T_F_uncr(:,j) = [T_F1_uncr';T_F2_uncr'];
    end
    for i = 1:size(obj.T_F,1)
        obj.T_F{i} = MEASvar(T_F(i,:),T_F_uncr(i,:),'C');
    end
    
    obj.pos_F = [22.7338*25.4e-3+0*y_F;y_F];

    group = 'Pressures_Inst';
    obj.m_dot_H = TDMSvar(tdms_struct,group,'Flow1'); 
    obj.m_dot_C = obj.m_dot_H;
    obj.P_H_in = TDMSvar(tdms_struct,group,'GP4');    
    obj.DP_H = TDMSvar(tdms_struct,group,'DP4');
    obj.P_C_in = TDMSvar(tdms_struct,group,'GP3');
    obj.DP_C = TDMSvar(tdms_struct,group,'DP3');


    function obj = TDMSvar(tdms_struct,group,variable)
        obj = struct('unit','-','data',1,'uncr',1,'mean',1,'merr',1);
        obj.data = tdms_struct.(group).(variable).data;
        obj.uncr = tdms_struct.([group,'_Error']).(variable).data;
        obj.unit = tdms_struct.(group).(variable).Props.unit_string;
        obj.mean = mean(obj.data);
        merr = sum(obj.uncr.^2)+sum(obj.data.^2)-length(obj.data)*obj.mean^2;
        obj.merr = sqrt( abs(merr)/length(obj.data) );
    end

    function obj = MEASvar(data,uncr,unit)
        obj = struct('unit','-','data',1,'uncr',1,'mean',1,'merr',1);
        obj.data = data;
        obj.uncr = uncr;
        obj.unit = unit;
        obj.mean = mean(obj.data);
        merr = sum(obj.uncr.^2)+sum(obj.data.^2)-length(obj.data)*obj.mean^2;
        obj.merr = sqrt( abs(merr)/length(obj.data) );
    end

end