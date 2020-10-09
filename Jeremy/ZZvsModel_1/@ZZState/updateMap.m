function hfig = updateMap(obj,hfig,seconds)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here
    
    nF = length(obj.Fiber.subsecInd);
    nP = length(obj.Fiber.subsecInd{1});
    TMap = zeros(nF,nP);

    
    if nargin == 2
    
        c_headers = [obj.T_H_in.mean;obj.T_H_out.mean;obj.T_C_in.mean;obj.T_C_out.mean];
        c_body = mean(c_headers);
        meanT = mean(obj.Fiber.Temperature);
        for i = 1:nF
            TMap(i,:) = meanT(1,obj.Fiber.subsecInd{i});
        end
        title_string = 'state at mean values';
    
    else
       
        dsec = seconds;
        fsec = dsec+obj.sec_diff;
        if dsec < obj.seconds(1)
            d_ind = 0;
        else
            d_ind = find((dsec<obj.seconds),1)-1;
        end
        
        if fsec < obj.Fiber.elapsedSeconds(1)
            f_ind = 0;
        else
            f_ind = find((fsec<obj.Fiber.elapsedSeconds),1)-1;
        end
        
        fnm = 1; dnm = 1;
        if isempty(d_ind) || (d_ind == 0)
            d_ind = 1; 
            dnm = NaN;
        end
        if isempty(f_ind) || (f_ind == 0)
            f_ind = 1; 
            fnm = NaN;
        end
        
        c_headers = dnm*[obj.T_H_in.data(d_ind);obj.T_H_out.data(d_ind);obj.T_C_in.data(d_ind);obj.T_C_out.data(d_ind)];
        c_body = mean(c_headers);
        matrT = obj.Fiber.Temperature(f_ind,:);
        for i = 1:nF
            TMap(i,:) = fnm*matrT(1,obj.Fiber.subsecInd{i});
        end
        title_string = ['state at ',num2str(seconds),' sec'];
        
    end
    
    figure(hfig)
    hfig.Children(2).Children(3).ZData = TMap;
    hfig.Children(2).Children(2).CData = c_headers;
    hfig.Children(2).Children(1).CData = c_body;
    
    title(title_string)
    
end