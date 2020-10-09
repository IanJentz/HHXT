function hfig = initMap(obj)
    %METHOD1 Summary of this method goes here
    %   Detailed explanation goes here
    
    c_headers = [obj.T_H_in.mean;obj.T_H_out.mean;obj.T_C_in.mean;obj.T_C_out.mean];
    c_body = mean(c_headers);
    
    nF = length(obj.Fiber.subsecInd);
    nP = length(obj.Fiber.subsecInd{1});
    
    meanT = mean(obj.Fiber.Temperature);
    TMap = zeros(nF,nP);
    for i = 1:nF
        TMap(i,:) = meanT(1,obj.Fiber.subsecInd{i});
    end
    
    title_string = 'state at mean values';
    
    T_max = max([ ...%max(max(obj.Fiber.Temperature)),...
        max(obj.T_H_in.data),max(obj.T_H_out.data),...
        max(obj.T_C_in.data),max(obj.T_C_out.data) ]);
    T_min = min([ ...%min(min(obj.Fiber.Temperature)),...
        min(obj.T_H_in.data),min(obj.T_H_out.data),...
        min(obj.T_C_in.data),min(obj.T_C_out.data) ]);
    
    hfig = figure();
    contourf(obj.pos_Map.x,obj.pos_Map.y,TMap)
%     levels = hfig.Children(1).Children(1).LevelList;
%     levels = [T_min,levels,T_max];
%     contourf(obj.pos_Map.x,obj.pos_Map.y,TMap,levels)
    caxis([T_min,T_max]);
    hfig.Children(1).Children(1).LineStyle = 'none';
    hold on
    patch(obj.x_headers,obj.y_headers,c_headers)
    patch(obj.x_body,obj.y_body,c_body,'FaceColor','none')
    colorbar    
    hold off
    
    hfig.OuterPosition = [30 550 1600 500];
    axis equal
    title(title_string,'interpreter','latex')
    hfig.OuterPosition = [30 550 1500 500];
    
end
