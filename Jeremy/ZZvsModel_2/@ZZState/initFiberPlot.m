function hfig = initFiberPlot(obj)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here
    
    nF = length(obj.Fiber.subsecInd);
    nP = length(obj.Fiber.subsecInd{1});
    
    meanT = mean(obj.Fiber.Temperature);
    TMap = zeros(nF,nP);
    for i = 1:nF
        TMap(i,:) = meanT(1,obj.Fiber.subsecInd{i});
    end
    title_string = 'state at mean values';
    
    hfig = figure();
        
    for i = 1:nF
        subplot(nF,1,i)
        hold on
        plot(obj.pos_Map.x,meanT(1,obj.Fiber.subsecInd{i}));
        plot([obj.pos_TF(1,i),obj.pos_TF(1,i+nF)],...
            [obj.T_F{i}.mean,obj.T_F{i+nF}.mean],'o');
        hold off
        if i == 1
            title(title_string,'interpreter','latex')
        end
    end
    hfig.OuterPosition = [500 50 600 1000]; 

    
end
