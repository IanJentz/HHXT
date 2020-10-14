function fig1 = PlotResults(filename,dir_plot,breakpoints)
% PlotResults - Plot and save various figures.  Creates plots for the
% convergence of the model, the temperature d.o.f's of the solution, the
% heating in W of the fluids, the fluid temperature distribution, the fluid
% pressure drop and velocity vectors, and the solid body temperature
% distribution.  These are saved as .fig, .pdf, and .png in both a large
% and small size.

if nargin == 1
    breakpoints = false;
    dir_plot = cd;
end
if nargin == 2
    breakpoints = false;
end
warning off
mkdir(dir_plot);
warning on

load(filename,'results','model')
T_C_in = 50;
T_H_in = 650;
T_C_in = results.Tables{3}.Temperature(1);
T_H_in = results.Tables{4}.Temperature(1);


fig1 = figure();
%plot convergence
cind = results.ConvStepData(:,2) == 1;
cdata = results.ConvStepData(cind,:);
steps = 1:sum(cind);

i = 4:sum(cind);
dcdata = zeros(size(cdata)); ddcdata = zeros(size(cdata));
dcdata(i,:) = cdata(i,:) - cdata(i-1,:);
ddcdata(i,:) = cdata(i-1,:) - 2*cdata(i,:) + cdata(i-1,:);
dcmean = mean(dcdata(:,3:9),2);
ddcmean = mean(ddcdata(:,3:9),2);
dcmean(1:(i(1)-1)) = dcmean(i(1)); ddcmean(1:(i(1)-1)) = ddcmean(i(1));
% dcmean(3:4) = 0; ddcmean(3:4) = 0;

idplot = 1:steps(end);

plot(steps,cdata(:,3:9)); set(gca,'YScale','log')

xlabel('number of full steps','Interpreter','latex')
ylabel('convergence value','Interpreter','latex')
lgd = legend('Residual','Abs. $\Delta$: Body $T$','Abs. $\Delta$: Fluid $P$','Abs. $\Delta$: Fluid $T$',...
    'Rel. $\Delta$: Body $T$','Rel. $\Delta$: Fluid $P$','Rel. $\Delta$: Fluid $T$');
 lgd.Interpreter = 'latex';
 lgd.NumColumns = 1;

fig1.Units = 'inches'; fig1.RendererMode = 'manual'; fig1.Renderer = 'painters';
fig1.OuterPosition = [2,3,3.34,5.45];
saveas(fig1,[dir_plot,'/Convergence.fig'],'fig');
saveas(fig1,[dir_plot,'/Convergence_334x545.png'],'png');  
saveas(fig1,[dir_plot,'/Convergence_334x545.pdf'],'pdf');
    
if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% create plotting object
Plot1 = hhxt.PlotHHXT(filename);
fig1 = gcf; fig1.Units = 'inches'; fig1.RendererMode = 'manual'; fig1.Renderer = 'painters';

%update ploting object:

% plot temperatures
Plot1.updatePlot('Temperatures','Wireframe','on','TempLimits',[T_C_in,T_H_in]);
fig1.OuterPosition = [2,3,6.85,7.5];
fig1.Children(9).Title.String = replace(fig1.Children(9).Title.String,'Body','Recuperator Body');
fig1.Children(7).Title.String = replace(fig1.Children(7).Title.String,'Fluid 1','Cold CO$_2$');
fig1.Children(5).Title.String = replace(fig1.Children(5).Title.String,'Fluid 2','Hot CO$_2$');
figure(fig1)
saveas(fig1,[dir_plot,'/Temps.fig'],'fig');
saveas(fig1,[dir_plot,'/Temps_685x750.png'],'png'); 
saveas(fig1,[dir_plot,'/Temps_685x750.pdf'],'pdf');
fig1.OuterPosition = [2,3,3.34,4];
saveas(fig1,[dir_plot,'/Temps_334x400.png'],'png'); 
saveas(fig1,[dir_plot,'/Temps_334x400.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot the heat going into each stream
Plot1.updatePlot('FluidHeating','Wireframe','on');
Cscale = 0.5*max(abs(min(min(results.StreamVolHeating))),abs(max(max(results.StreamVolHeating))));
Plot1.updatePlot('FluidHeating','Wireframe','on','CLimits',Cscale*[-1,1]); % fix the limits of the color bar
Cfact = model.Mesh.MaxElementSize/0.01;
Plot1.updatePlot('FluidHeating','Wireframe','on','CLimits',Cfact*max(max(abs(results.StreamVolHeating)))*[-1,1]);
fig1.OuterPosition = [2,3,6.85,5.4];
fig1.Children(6).Title.String = replace(fig1.Children(6).Title.String,'Fluid 1','Cold CO$_2$');
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 2','Hot CO$_2$');
figure(fig1)
saveas(fig1,[dir_plot,'/Heating.fig'],'fig');
saveas(fig1,[dir_plot,'/Heating_685x540.png'],'png');  
saveas(fig1,[dir_plot,'/Heating_685x540.pdf'],'pdf');
fig1.OuterPosition = [2,3,3.34,3];
saveas(fig1,[dir_plot,'/Heating_334x300.png'],'png');  
saveas(fig1,[dir_plot,'/Heating_334x300.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot fluid temperatures
Plot1.updatePlot('FluidTemperatures','Wireframe','on','TempLimits',[T_C_in,T_H_in]);
fig1.OuterPosition = [2,3,6.85,5];
fig1.Children(6).Title.String = replace(fig1.Children(6).Title.String,'Fluid 1','Cold CO$_2$');
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 2','Hot CO$_2$');
figure(fig1)
saveas(fig1,[dir_plot,'/FluidTemps.fig'],'fig');
saveas(fig1,[dir_plot,'/FluidTemps_685x500.png'],'png'); 
saveas(fig1,[dir_plot,'/FluidTemps_685x500.pdf'],'pdf');
fig1.OuterPosition = [2,3,3.34,3];
saveas(fig1,[dir_plot,'/FluidTemps_334x300.png'],'png'); 
saveas(fig1,[dir_plot,'/FluidTemps_334x300.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

% plot pressure and velocity vectors
Plot1.updatePlot('Press+ChannelVelocityVector','Wireframe','on','pdeplotArgs',{'AutoScaleFactor',2,'Color','magenta'});
fig1.OuterPosition = [2,3,6.85,5.6];
fig1.Children(6).Title.String = replace(fig1.Children(6).Title.String,'Fluid 1','Cold CO$_2$');
fig1.Children(4).Title.String = replace(fig1.Children(4).Title.String,'Fluid 2','Hot CO$_2$');
figure(fig1)
saveas(fig1,[dir_plot,'/VelocityVectors.fig'],'fig');
saveas(fig1,[dir_plot,'/VelocityVectors_685x560.png'],'png'); 
saveas(fig1,[dir_plot,'/VelocityVectors_685x560.pdf'],'pdf');
fig1.OuterPosition = [2,3,3.34,3.6];
saveas(fig1,[dir_plot,'/VelocityVectors_334x360.png'],'png'); 
saveas(fig1,[dir_plot,'/VelocityVectors_334x360.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

%plot the solid body temperature
Plot1.updatePlot('BodyTemperature','Mesh','on','TempLimits',[T_C_in,T_H_in]);
fig1.OuterPosition = [2,3,6.85,3.3];
fig1.Children(3).Title.String = replace(fig1.Children(3).Title.String,'Body','Recuperator Body');
figure(fig1)
saveas(fig1,[dir_plot,'/BodyTemp.fig'],'fig');
saveas(fig1,[dir_plot,'/BodyTemp_685x330.png'],'png'); 
saveas(fig1,[dir_plot,'/BodyTemp_685x330.pdf'],'pdf');
fig1.OuterPosition = [2,3,3.34,2.3];
saveas(fig1,[dir_plot,'/BodyTemp_334x230.png'],'png'); 
saveas(fig1,[dir_plot,'/BodyTemp_334x230.pdf'],'pdf');


if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

end