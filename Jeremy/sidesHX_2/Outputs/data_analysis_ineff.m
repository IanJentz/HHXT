%% Side headers data analysis
%load the results into the workspace
results = readtable('finaloutput.csv');
close all
%since there are three DOF, results X will be plotted for one variable A
%fixed at a value close to the original case, and the second variable B
%varying for different values of the third variable C
%e.g. plot ineff (X) as a function of alpha (A) for 3 values of C_r (C)
%with NTU (B) fixed at 5
%values for which variable C will plotted
C_r_index = [0.1 0.3 0.5 0.7 1];
NTU_index = [1.5849 6.3096 25.119 63.096];
alpha_index = [0 0.2 0.5 0.7 1];

%% plotting the effect of Area ratio

%fixed NTU, various C_r
figure;
for i=1:length(C_r_index)
    cond = results.C_r == C_r_index(i) & results.NTU == 6.3096;
    ineff = results.ineff(cond);
    alpha = results.alpha(cond);
    plot(alpha,ineff);
    hold on
end
title('Effect of Area Ratio on Ineffectivness for Various C$_{r}$','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(C_r_index(:)))
title(legend,'C_r')
axis([0 1 -0.1 0.5])
saveas(figure(1),['figures','/','alpha for various C_r','.png'],'png');  

%fixed C_r, various NTU
figure;
for i=1:length(NTU_index)
    cond = results.NTU == NTU_index(i) & results.C_r == 1;
    ineff = results.ineff(cond);
    alpha = results.alpha(cond);
    plot(alpha,ineff);
    hold on
end
title('Effect of Area ratio on Ineffectivness for Various NTU','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)))
title(legend,'NTU')
axis([0 1 -0.1 0.5])
saveas(figure(2),['figures','/','alpha for various NTU','.png'],'png');  

%% plotting the effect of NTU

%fixed alpha, various C_r
figure;
for i=1:length(C_r_index)
    cond = results.C_r == C_r_index(i) & results.alpha == 0.2;
    ineff = results.ineff(cond);
    NTU = results.NTU(cond);
    data = sortrows([NTU ineff],1);
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of NTU on Ineffectivness for Various C$_{r}$','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(C_r_index(:)))
title(legend,'C_r')
set(gca, 'XScale', 'log')
axis([1 100 -0.1 0.5])
saveas(figure(3),['figures','/','NTU for various C_r','.png'],'png');  

%fixed C_r, various alpha
figure;
for i=1:length(alpha_index)
    cond = results.alpha == alpha_index(i) & results.C_r == 1;
    ineff = results.ineff(cond);
    NTU = results.NTU(cond);
    data = sortrows([NTU ineff],1);
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of NTU on Ineffectivness for Various $\alpha$ and $C_r$=1','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)))
title(legend,'alpha')
set(gca, 'XScale', 'log')
axis([1 100 -0.1 0.5])
saveas(figure(4),['figures','/','NTU for various alpha C_r 1','.png'],'png');  


% 
% figure;
% for i=1:length(alpha_index)
%     cond = results.alpha == alpha_index(i) & results.C_r == 0.7;
%     ineff = results.ineff(cond);
%     NTU = results.NTU(cond);
%     data = sortrows([NTU ineff],1);
%     plot(data(:,1),data(:,2));
%     hold on
% end
% title('Effect of NTU on Ineffectivness for Various $\alpha$ and $C_r$=0.7','interpreter','latex')
% xlabel('NTU $\left[ - \right]$','interpreter','latex')
% ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
% grid on
% legend(num2str(alpha_index(:)))
% title(legend,'alpha')
% set(gca, 'XScale', 'log')
% axis([1 100 -0.1 0.5])
% saveas(figure(1),['figures','/','NTU for various alpha C_r 0.7','.png'],'png');  
% 


%% plotting the effect of C_r

%fixed alpha, various NTU
figure;
for i=1:length(NTU_index)
    cond = results.NTU == NTU_index(i) & results.alpha == 0.2;
    ineff = results.ineff(cond);
    C_r = results.C_r(cond);
    plot(C_r,ineff);
    hold on
end
title('Effect of Capacitance Ratio on Ineffectivness for Various NTU','interpreter','latex')
xlabel('Capacitance Ratio C$_{r}$ $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)))
title(legend,'NTU')
axis([0 1 -0.1 0.5])
saveas(figure(5),['figures','/','C_r for various NTU','.png'],'png');  

%fixed NTU, various alpha
figure;
for i=1:length(alpha_index)
    cond = results.alpha == alpha_index(i) & results.NTU == 6.3096;
    ineff = results.ineff(cond);
    C_r = results.C_r(cond);
    plot(C_r,ineff);
    hold on
end
title('Effect of Capacitance Ratio on Ineffectivness for Various $\alpha$','interpreter','latex')
xlabel('Capacitance Ratio C$_{r}$ $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)))
title(legend,'alpha')
axis([0 1 -0.1 0.5])
saveas(figure(6),['figures','/','C_r for various alpha','.png'],'png');  
