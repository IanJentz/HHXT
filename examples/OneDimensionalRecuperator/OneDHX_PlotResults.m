figdir = 'figures'; %figure directory to write into
mkdir(figdir);
pathdelim = '\';

%% plot results
%interpolate the result from the centerline
T_center = InterpolateLines(results,L);

%plot and save temperature distribution
hfig = figure();
hold on 
plot(T_center.x,T_center.T_H,'r');
plot(T_center.x,T_center.T_C,'b');
plot(T_center.x,T_center.T_s,'g');
hold off
xlabel('position $\left[ m \right]$','interpreter','latex')
ylabel('temperature $\left[ C \right]$','interpreter','latex')
% lgd = legend('EES : $T_H$','EES : $T_C$');
lgd = legend('$T_H$','$T_C$','$T_S$');
lgd.Interpreter = 'latex';
filename = 'temperatures';
saveas(hfig,[figdir,pathdelim,filename,'.fig'],'fig');  
saveas(hfig,[figdir,pathdelim,filename,'.png'],'png');  
saveas(hfig,[figdir,pathdelim,filename,'.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint
end

%plot and save pressure distribution
hfig = figure();
hold on 
plot(T_center.x,T_center.P_H,'r');
plot(T_center.x,T_center.P_C,'b');
hold off
xlabel('position $\left[ m \right]$','interpreter','latex')
ylabel('pressure $\left[ Pa \right]$','interpreter','latex')
lgd = legend('$P_H$','$P_C$');
lgd.Interpreter = 'latex';
filename = 'pressures';
saveas(hfig,[figdir,pathdelim,filename,'.fig'],'fig');  
saveas(hfig,[figdir,pathdelim,filename,'.png'],'png');  
saveas(hfig,[figdir,pathdelim,filename,'.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint
end

%plot and save velocity distribution
hfig = figure();
hold on 
plot(T_center.x,T_center.v_H,'r');
plot(T_center.x,T_center.v_C,'b');
hold off
xlabel('position $\left[ m \right]$','interpreter','latex')
ylabel('channel velocity $\left[ \frac{m}{s} \right]$','interpreter','latex')
lgd = legend('$\mathbf{v}_H$','$\mathbf{v}_C$');
lgd.Interpreter = 'latex';
filename = 'velocities';
saveas(hfig,[figdir,pathdelim,filename,'.fig'],'fig');  
saveas(hfig,[figdir,pathdelim,filename,'.png'],'png');  
saveas(hfig,[figdir,pathdelim,filename,'.pdf'],'pdf');

if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint
end

