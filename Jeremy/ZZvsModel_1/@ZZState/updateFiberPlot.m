function hfig = updateFiberPlot(obj,hfig,seconds)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here
    
    nF = length(obj.Fiber.subsecInd);
    nP = length(obj.Fiber.subsecInd{1});
    TMap = zeros(nF,nP);

    
    if nargin == 2
    
        meanT = mean(obj.Fiber.Temperature);
        TMap = zeros(nF,nP);
        for i = 1:nF
            TMap(i,:) = meanT(1,obj.Fiber.subsecInd{i});
        end
        title_string = 'state at mean values';
        
        figure(hfig)
        for i = 1:nF
            hfig.Children(nF-i+1).Children(2).YData = meanT(1,obj.Fiber.subsecInd{i});
            hfig.Children(nF-i+1).Children(1).YData = [obj.T_F{i}.mean,obj.T_F{i+nF}.mean];        

            if i == 1
               hfig.Children(nF-i+1).Title.String = title_string;
            end
        end
    
    else
       
        dsec = seconds;
        fsec = dsec+obj.sec_diff;
        d_ind = find((dsec<obj.seconds),1)-1;
        f_ind = find((fsec<obj.Fiber.elapsedSeconds),1)-1;
        dnm = 1; fnm = 1;
        if isempty(d_ind) || (d_ind == 0)
            d_ind = 1; 
            dnm = NaN;
        end
        if isempty(f_ind) || (f_ind == 0)
            f_ind = 1; 
            fnm = NaN;
        end
        
        matrT = fnm*obj.Fiber.Temperature(f_ind,:);
        title_string = ['state at ',num2str(seconds),' sec'];
        
        figure(hfig)
        for i = 1:nF
            hfig.Children(nF-i+1).Children(2).YData = matrT(1,obj.Fiber.subsecInd{i});
try
            hfig.Children(nF-i+1).Children(1).YData = dnm*[obj.T_F{i}.data(d_ind),obj.T_F{i+nF}.data(d_ind)];        
catch
    keyboard
end
            if i == 1
               hfig.Children(nF-i+1).Title.String = title_string;
            end
        end
        
    end
    
    
end