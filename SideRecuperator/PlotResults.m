function fig1 = PlotResults(filename,dir_plot,breakpoints,R_V_cond,strm)

sfigfile = false;
spdffile = false;

if nargin == 1
    breakpoints = false;
    dir_plot = cd;
end
if nargin == 2
    breakpoints = false;
end
warning off
mkdir([dir_plot,'/figs']);
warning on

load(filename,'results','model')
T_C_in = results.Tables{3}.Temperature(1);
T_H_in = results.Tables{4}.Temperature(1);

fNames = {'Cold CO$_2$','Hot CO$_2$'};

Tlims = [T_C_in,T_H_in];
Tlims = [min(Tlims),max(Tlims)];

fig1 = figure();
warning off
mkdir('figs');
warning on
%plot convergence
% cind = model.Runtime.ConvergenceData(:,2) == 1;
% cdata = model.Runtime.ConvergenceData(cind,:);
cind = results.ConvStepData(:,2) == 1;
cdata = results.ConvStepData(cind,:);
steps = 1:sum(cind);

if sum(cind)>=2
i = 2:sum(cind);
dcdata = zeros(size(cdata)); ddcdata = zeros(size(cdata));
dcdata(i,:) = cdata(i,:) - cdata(i-1,:);
ddcdata(i,:) = cdata(i-1,:) - 2*cdata(i,:) + cdata(i-1,:);
dcmean = mean(dcdata(:,3:9),2);
ddcmean = mean(ddcdata(:,3:9),2);
dcmean(1:(i(1)-1)) = dcmean(i(1)); ddcmean(1:(i(1)-1)) = ddcmean(i(1));
% dcmean(3:4) = 0; ddcmean(3:4) = 0;
end

idplot = 1:steps(end);

% subplot(3,1,1)
p1 = plot(steps,cdata(:,3:9)); set(gca,'YScale','log')
p1(1).Color = '#000000';
p1(2).Color = '#FF00FF';
p1(3).Color = '#0000FF';
p1(4).Color = '#FF0000';
p1(5).Color = '#FF00FF';
p1(6).Color = '#0000FF';
p1(7).Color = '#FF0000';
p1(1).LineStyle = '-';
p1(2).LineStyle = ':';
p1(3).LineStyle = ':';
p1(4).LineStyle = ':';
p1(5).LineStyle = '--';
p1(6).LineStyle = '--';
p1(7).LineStyle = '--';

% subplot(3,1,2)
% plot(idplot,dcmean(idplot));
% subplot(3,1,3)
% plot(idplot,ddcmean(idplot));

AbsTol_P = 10*model.SolverOptions.AbsoluteTolerance;
AbsTol_T = model.SolverOptions.AbsoluteTolerance;
RelTol_P = model.SolverOptions.RelativeTolerance;
RelTol_T = model.SolverOptions.RelativeTolerance;
ResTol = model.SolverOptions.ResidualTolerance;
Tols = [ResTol,AbsTol_T,AbsTol_P,AbsTol_T,RelTol_T,RelTol_P,RelTol_T];
lgditems = {'Residual','Abs. $\Delta$: Body $T$','Abs. $\Delta$: Fluid $P$','Abs. $\Delta$: Fluid $T$',...
    'Rel. $\Delta$: Body $T$','Rel. $\Delta$: Fluid $P$','Rel. $\Delta$: Fluid $T$'};
for k = 1:length(lgditems)
   lgditems{k} = [lgditems{k}, ', $\leq $',num2str(Tols(k),'%5.0e')];
end

xlabel('number of full steps','Interpreter','latex')
ylabel('convergence value','Interpreter','latex')
lgd = legend(lgditems,'Location','eastoutside');
 lgd.Interpreter = 'latex';
 lgd.NumColumns = 1;
 
fig1.Children(2).YGrid = 'on';
fig1.Children(2).YMinorGrid = 'on';

fig1.Units = 'inches'; fig1.RendererMode = 'manual'; fig1.Renderer = 'painters';
fig1.OuterPosition = [2,3,6.85,4.4];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
if sfigfile == true; saveas(fig1,[dir_plot,'/figs/Convergence.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/Convergence_2col.png'],'png');  
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Convergence_2col.pdf'],'pdf'); end
    
if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

close(fig1)

% create plotting object
Plot1 = hhxt.PlotHHXT(filename);
fig1 = gcf; fig1.Units = 'inches'; fig1.RendererMode = 'manual'; fig1.Renderer = 'painters';

%update ploting object:

% plot temperatures
Plot1.updatePlot('Temperatures','Wireframe','on','TempLimits',Tlims);
fig1.OuterPosition = [2,3,6.85,7.5];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(9).Title.String = replace(fig1.Children(9).Title.String,'Body','Recuperator Body');
fig1.Children(7).Title.String = replace(fig1.Children(7).Title.String,'Fluid 1',fNames{strm(1)});
fig1.Children(5).Title.String = replace(fig1.Children(5).Title.String,'Fluid 2',fNames{strm(2)});
figure(fig1)
if sfigfile == true; saveas(fig1,[dir_plot,'/figs/Temps.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/Temps_2col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Temps_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,4];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/Temps_1col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Temps_1col.pdf'],'pdf'); end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot the DT accross stream
Plot1.updatePlot('FluidDTs','Wireframe','on','R_V_cond',R_V_cond);
fig1.OuterPosition = [2,3,6.85,5.4];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
figure(fig1)
if sfigfile == true; saveas(fig1,[dir_plot,'/figs/Heating.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/StreamDTs_2col.png'],'png');  
if spdffile == true; saveas(fig1,[dir_plot,'/figs/StreamDTs_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,3];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/StreamDTs_1col.png'],'png');  
if spdffile == true; saveas(fig1,[dir_plot,'/figs/StreamDTs_1col.pdf'],'pdf'); end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot the heat going into each stream
Plot1.updatePlot('FluidHeating','Wireframe','on');
Cscale = max(abs(min(min(results.StreamVolHeating))),abs(max(max(results.StreamVolHeating))));
Plot1.updatePlot('FluidHeating','Wireframe','on','CLimits',Cscale*[-1,1]); % fix the limits of the color bar
Cfact = model.Mesh.MaxElementSize/0.01;
% Plot1.updatePlot('FluidHeating','Wireframe','on','CLimits',Cfact*max(max(abs(results.StreamVolHeating)))*[-1,1]);
fig1.OuterPosition = [2,3,6.85,5.4];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(6).Title.String = replace(fig1.Children(6).Title.String,'Fluid 1',fNames{strm(1)});
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 2',fNames{strm(2)});
figure(fig1)
if sfigfile == true; saveas(fig1,[dir_plot,'/figs/Heating.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/Heating_2col.png'],'png');  
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Heating_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,3];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/Heating_1col.png'],'png');  
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Heating_1col.pdf'],'pdf'); end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot fluid temperatures
Plot1.updatePlot('FluidTemperatures','Wireframe','on','TempLimits',Tlims);
fig1.OuterPosition = [2,3,6.85,5];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(6).Title.String = replace(fig1.Children(6).Title.String,'Fluid 1',fNames{strm(1)});
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 2',fNames{strm(2)});
figure(fig1)
if sfigfile == true; saveas(fig1,[dir_plot,'/figs/FluidTemps.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/FluidTemps_2col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/FluidTemps_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,3];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/FluidTemps_1col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/FluidTemps_1col.pdf'],'pdf'); end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot pressure and velocity vectors
Plot1.updatePlot('Press+ChannelVelocityVector','Wireframe','on','Elemental','pdeplotArgs',{'AutoScaleFactor',2,'Color','magenta'});
fig1.OuterPosition = [2,3,6.85,5.6];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(6).Title.String = replace(fig1.Children(6).Title.String,'Fluid 1',fNames{strm(1)});
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 2',fNames{strm(2)});
figure(fig1)
if sfigfile == true;  saveas(fig1,[dir_plot,'/figs/VelocityVectors.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/VelocityVectors_2col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/VelocityVectors_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,3.6];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/VelocityVectors_1col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/VelocityVectors_1col.pdf'],'pdf');end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

%plot the solid body temperature
Plot1.updatePlot('BodyTemperature','Mesh','on','TempLimits',Tlims);
fig1.OuterPosition = [2,3,6.85,3.3];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(3).Title.String = replace(fig1.Children(3).Title.String,'Body','Recuperator Body');
figure(fig1)
if sfigfile == true;  saveas(fig1,[dir_plot,'/figs/BodyTemp.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/BodyTemp_2col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/BodyTemp_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,2.3];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/BodyTemp_1col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/BodyTemp_1col.pdf'],'pdf');end


if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot pressure and velocity vectors
Cscale = [min(min(results.ElemMassFlux(results.ElemMassFlux~=0))),max(max(results.ElemMassFlux))];
Plot1.updatePlot('ChannelMassFluxMagnitude','Wireframe','on','Elemental','CLimits',Cscale);
fig1.OuterPosition = [2,3,6.85,5.6];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 1',fNames{strm(1)});
fig1.Children(2).Title.String = replace(fig1.Children(2).Title.String,'Fluid 2',fNames{strm(2)});
figure(fig1)
if sfigfile == true;  saveas(fig1,[dir_plot,'/figs/MassFlux.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/MassFlux_2col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/MassFlux_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,3.6];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/MassFlux_1col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/MassFlux_1col.pdf'],'pdf');end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot pressure and velocity vectors
Cscale = [min(min(results.ElemReynolds(results.ElemReynolds~=0))),max(max(results.ElemReynolds))];
Plot1.updatePlot('ReynoldsMagnitude','Wireframe','on','Elemental','CLimits',Cscale);
fig1.OuterPosition = [2,3,6.85,5.6];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 1',fNames{strm(1)});
fig1.Children(2).Title.String = replace(fig1.Children(2).Title.String,'Fluid 2',fNames{strm(2)});
figure(fig1)
if sfigfile == true;  saveas(fig1,[dir_plot,'/figs/Reynolds.fig'],'fig'); end
saveas(fig1,[dir_plot,'/figs/Reynolds_2col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Reynolds_2col.pdf'],'pdf'); end
fig1.OuterPosition = [2,3,3.34,3.6];
fig1.PaperSize = fig1.OuterPosition([end-1,end]);
saveas(fig1,[dir_plot,'/figs/Reynolds_1col.png'],'png'); 
if spdffile == true; saveas(fig1,[dir_plot,'/figs/Reynolds_1col.pdf'],'pdf');end

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

close(fig1)

end