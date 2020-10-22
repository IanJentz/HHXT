function obj = calcJ(obj,date_start,elapsed_sec)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here

    testing = obj.debug;

    % determine how much time to take
    if nargin == 1
        sub_ind = 1:obj.ntimes;
    else
       sec_start = elapsed(date_start-obj.time_start);
       i_s = find( obj.seconds >= sec_start ,1);
       if nargin == 2
           i_e = length(obj.data);
       else
           i_e = find( obj.seconds <= sec_start+elapsed_sec ,1);
       end
       sub_ind = i_s:1:i_e;
    end

    %write subset mean and mean error to file
    mdot = mean_err(obj.m_dot_H,sub_ind);
    GPH = mean_err(obj.P_H_in,sub_ind);
    DPH = mean_err(obj.DP_H,sub_ind);
    GPC = mean_err(obj.P_C_in,sub_ind);
    DPC = mean_err(obj.DP_C,sub_ind);
    TCHin = mean_err(obj.T_H_in,sub_ind);
    TCHout = mean_err(obj.T_H_out,sub_ind);
    TCCin = mean_err(obj.T_C_in,sub_ind);
    TCCout = mean_err(obj.T_C_out,sub_ind);
    
    file = [obj.wdir,'\AveData.dat'];
    fid = fopen(file, 'w')  ;  
    fprintf(fid,'Flow_exp\t%8f\n',mdot.mean);   
    fprintf(fid,'GPH_exp\t%8f\n',GPH.mean); 
    fprintf(fid,'DPH_exp\t%8f\n',DPH.mean); 
    fprintf(fid,'GPC_exp\t%8f\n',GPC.mean); 
    fprintf(fid,'DPC_exp\t%8f\n',DPC.mean); 
    fprintf(fid,'TCHin_exp\t%8f\n',TCHin.mean); 
    fprintf(fid,'TCHout_exp\t%8f\n',TCHout.mean); 
    fprintf(fid,'TCCin_exp\t%8f\n',TCCin.mean); 
    fprintf(fid,'TCCout_exp\t%8f\n',TCCout.mean);
    fclose(fid) ;
    file = [obj.wdir,'\AveError.dat'];
    fid = fopen(file, 'w')  ;  
    fprintf(fid,'Flow_exp_err\t%8f\n',mdot.merr);   
    fprintf(fid,'GPH_exp_err\t%8f\n',GPH.merr); 
    fprintf(fid,'DPH_exp_err\t%8f\n',DPH.merr); 
    fprintf(fid,'GPC_exp_err\t%8f\n',GPC.merr); 
    fprintf(fid,'DPC_exp_err\t%8f\n',DPC.merr); 
    fprintf(fid,'TCHin_exp_err\t%8f\n',TCHin.merr); 
    fprintf(fid,'TCHout_exp_err\t%8f\n',TCHout.merr); 
    fprintf(fid,'TCCin_exp_err\t%8f\n',TCCin.merr); 
    fprintf(fid,'TCCout_exp_err\t%8f\n',TCCout.merr);
    fclose(fid) ;

    %run EES code to calculate J
    file = [obj.wdir,'\Calculate_ZZ.bat'];
    if testing; file = [file(1:end-4),'_test',file(end-3:end)]; end
    [status,cmdout] = system(file,'-echo');

    %read in EES Solution
    file = [obj.wdir,'\calcSol.dat'];
    fid = fopen(file,'r');
    tline = fgetl(fid);
    fclose(fid);
    tline = strsplit(tline,'-');
    cols = str2double(tline{1});
    rows = str2double(tline{2});

    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [1, 103];
    opts.Delimiter = [" ", "="]; 
    opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
    opts.VariableTypes = ["string", "double", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, "VarName1", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["VarName1", "VarName3"], "EmptyFieldRule", "auto");
    tbl = readtable(file,opts);
    VarName1 = tbl.VarName1;
    VarName2 = tbl.VarName2;
    VarName3 = tbl.VarName3;

    argin = cell(1,2*length(VarName1)); argin2 = argin;
    ind = [];
    for i = 1:length(VarName1)
        if VarName1{i}(end)=='$'
            ind = [ind,2*i-1,2*i];
        else
            argin{2*i-1} = string(VarName1(i));
            argin{2*i} = VarName2(i);
            argin2{2*i-1} = string(VarName1(i));
            argin2{2*i} = string(VarName3(i));
        end
    end
    argin(ind) = []; argin2(ind) = []; 
    obj.EES_sol = struct(argin{:});
    obj.EES_sol_units = struct(argin2{:});

    %read in EES arrays
    file = [obj.wdir,'\calcArr.dat'];

    fid = fopen(file,'r');
    tline = fgetl(fid);
    fclose(fid);
    tline = strsplit(tline,'-');
    cols = str2double(tline{1});
    rows = str2double(tline{2});

    opts = delimitedTextImportOptions("NumVariables", 28);
    opts.DataLines = [2, rows+1];
    opts.Delimiter = ["\t", " "];
    opts.VariableNames = ["Var1", "VarName2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"];
    opts.SelectedVariableNames = ["VarName2", "VarName3"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Var1", "VarName2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "VarName2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28"], "EmptyFieldRule", "auto");
    tbl = readtable(file,opts);
    VarName2 = tbl.VarName2;
    VarName3 = tbl.VarName3;

    argin = cell(1,2*length(VarName2));
    names = cell(1,length(VarName2)); units = names;
    for i = 1:length(VarName2)
            names{i} = string(VarName2{i}(1:end-3));
            argin{2*i-1} = names{i};
            units{i} = string(VarName3{i}(2:end-1));
            argin{2*i} = units{i};
    end
%             obj.EES_arr_units = struct(argin{:});

    opts = delimitedTextImportOptions("NumVariables", rows);
    opts.DataLines = [rows+2, Inf];
    opts.Delimiter = ["\t", " "];
    opts.VariableNames = string(names);
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    calcArr = readtable("W:\PCHE\AFandZZ\postProcessing\calcArr.dat", opts);
    calcArr.Properties.VariableUnits = string(units);
    obj.EES_arr = calcArr;

    %read in EES arrays
    file = [obj.wdir,'\calcErr.dat'];

    fid = fopen(file,'r');
    tline = fgetl(fid);
    fclose(fid);
    tline = strsplit(tline,'-');
    cols = str2double(tline{1});
    rows = str2double(tline{2});

    opts = delimitedTextImportOptions("NumVariables", 38);

    opts.DataLines = [2, rows+1];
    opts.Delimiter = ["\t", " "];
    opts.VariableNames = ["Var1", "VarName2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38"];
    opts.SelectedVariableNames = ["VarName2", "VarName3"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, ["Var1", "VarName2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "VarName2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38"], "EmptyFieldRule", "auto");
    tbl = readtable("W:\PCHE\AFandZZ\postProcessing\calcErr.dat", opts);
    VarName2 = tbl.VarName2;
    VarName3 = tbl.VarName3;

    argin = cell(1,2*length(VarName2));
    names = cell(1,length(VarName2)); units = names;
    for i = 1:length(VarName2)
            name = string(VarName2(i));
            if (contains(VarName2{i},'err'))
                name = char(name);
                name = string(name(1:end-4));
            end
            names{i} = name;
            argin{2*i-1} = name;
            units{i} = string(VarName3(i));
            argin{2*i} = units{i};
    end

    fid = fopen(file,'r');
    for i = 1:rows+2
    tline = fgetl(fid);
    end
    fclose(fid);
    pline = strsplit(tline,'\t');
    tline = strsplit(tline,{'\t','±'});
    tline = str2double(tline);

    j = 1;
    values = cell(1,length(tline)); 
    varnames = values;
    varunits = values;
    for i=1:length(pline)  
        varnames{j} = char(names{i});
        varunits{j} = char(units{i});
        values{j} = tline(j);
        j = j+1;
        if ismember('±',pline{i})
        varnames{j} =  [char(names{i}(1,1)),'_err'];
        varunits{j} = char(units{i});
        values{j} = tline(j);
        j = j+1; 
        end
    end
    [varnames,ia,~] = unique(varnames);
    varunits = varunits(ia);
    values = values(ia);
    tbl = table(values{:});
    for i = 1:length(varnames)
        tbl.Properties.VariableNames{i} = varnames{i};
        tbl.Properties.VariableUnits{i} = varunits{i};
    end

    obj.EES_err = tbl;


    function obj = mean_err(obj,sub_ind)
        obj.mean = mean(obj.data(sub_ind));
        merr = sum(obj.uncr(sub_ind).^2)+sum(obj.data(sub_ind).^2)-length(obj.data(sub_ind))*obj.mean^2;
        obj.merr = sqrt( abs(merr)/length(obj.data(sub_ind)) );
    end
end
        