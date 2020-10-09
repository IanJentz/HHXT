function [cmp_fiber,cmp_TC,fpnt,epnt] = ModelVSFibers(obj,results,dir_plot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if nargin < 3
    plotting = false;
    dir_plot = cd;
else
    plotting = true;
end

% for ZZ recuperator
L = 22.7338*25.4e-3; W = 4.6774*25.4e-3; H = 33e-3;
W_head = 3.3865*25.4e-3;
W_chan = 3.6777*25.4e-3;

% pull in fiber data
nF = length(obj.Fiber.subsecInd);
nP = length(obj.Fiber.subsecInd{1});

fiber_x = obj.pos_Map.x - 0.5*sum(obj.pos_Map.x([1,end]));
fiber_y = obj.pos_Map.y - 0.5*sum(obj.pos_Map.y([4,5]));

meanT = mean(obj.Fiber.Temperature);
T_model = zeros(nF,length(fiber_x));
T_fiber = T_model;
T_TC = zeros(nF,2);
T_TC_x = T_TC;
T_TC_model = T_TC;
for i = 1:nF
    T_TC(i,:) = [obj.T_F{i}.mean,obj.T_F{i+nF}.mean];
    T_TC_x(i,:) = [obj.pos_TF(1,i),obj.pos_TF(1,i+nF)]-0.5*sum(obj.pos_Map.x([1,end]));
    T_fiber(i,:) = meanT(1,obj.Fiber.subsecInd{i});
    T_model(i,:) = interpolateSolution(results,fiber_x,fiber_y(i)*ones(1,length(fiber_x)),1);
    T_TC_model(i,:) = interpolateSolution(results,T_TC_x(i,:),fiber_y(i)*ones(1,2),1);
end

dif_TC = T_TC_model-T_TC; ddif_TC = abs(dif_TC(:));
cmp_TC = [rms(dif_TC),std(dif_TC)];
dif_fiber = T_model-T_fiber; dif_fiber = abs(dif_fiber(:));
dif_fiber(isnan(dif_fiber)) = [];
cmp_fiber = [rms(dif_fiber),std(dif_fiber)];

inF = 4;
ix = floor(length(fiber_x)/2);
fpnt = [T_model(inF,ix),T_TC_model(inF,1)];
epnt = [T_fiber(inF,ix),T_TC(inF,1)];

if plotting == true
    warning off
    mkdir([dir_plot,'/figs']);
    warning on
    
    title_string = 'state at mean values';

    fig1 = figure();
        
    for i = 1:nF
        subplot(nF,1,i)
        hold on
        plot(fiber_x,T_fiber(i,:),'k');
        plot(fiber_x,T_model(i,:),'r')
        plot(T_TC_x(i,:),T_TC(i,:),'ok');
        hold off
        if i == 1
            title(title_string,'interpreter','latex')
        end
    end
    fig1.OuterPosition = [500 50 600 1000]; 
    
    fig1.Units = 'inches'; fig1.RendererMode = 'manual'; fig1.Renderer = 'painters';
    fig1.OuterPosition = [2,1,3.34,9.60];
    saveas(fig1,[dir_plot,'/figs/ModelComp_1.fig'],'fig');
    saveas(fig1,[dir_plot,'/figs/ModelComp_1_334x960.png'],'png');  
    saveas(fig1,[dir_plot,'/figs/ModelComp_1_334x960.pdf'],'pdf');
    
    
       
    
end


end

