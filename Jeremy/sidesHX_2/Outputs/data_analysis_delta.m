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
alpha_index = [0.2 0.5 0.7 1];

%% plotting the effect of Area ratio

%fixed NTU, various C_r
figure;
for i=1:length(C_r_index)
    cond = results.C_r == C_r_index(i) & results.NTU == 6.3096;
    eff_cf = results.eff(cond&results.alpha == 0);
    eff = results.eff(cond);
    DELTAeff = eff_cf-eff;
    alpha = results.alpha(cond);
    plot(alpha,DELTAeff);
    hold on
end
title('Effect of Area Ratio on Effectiveness for Various C$_{r}$ (NTU=6.3096)','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(C_r_index(:)),'Location','northwest')
title(legend,'C_r')
axis([0 1 -0.02 0.1])
saveas(figure(1),['figures_eff','/','alpha for various C_r','.png'],'png');  

%fixed C_r, various NTU
figure;
for i=1:length(NTU_index)
    cond = results.NTU == NTU_index(i) & results.C_r == 1;
    eff_cf = results.eff(cond&results.alpha == 0);
    eff = results.eff(cond);
    DELTAeff = eff_cf-eff;
    alpha = results.alpha(cond);
    plot(alpha,DELTAeff);
    hold on
end
title('Effect of Area ratio on Effectiveness for Various NTU ($C_r$=1)','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)),'Location','northwest')
title(legend,'NTU')
axis([0 1 -0.02 0.1])
saveas(figure(2),['figures_eff','/','alpha for various NTU','.png'],'png');  

%% plotting the effect of NTU

%fixed alpha, various C_r
figure;
for i=1:length(C_r_index)
    cond = results.C_r == C_r_index(i) & results.alpha == 0.2;
    cond_cf = results.C_r == C_r_index(i) & results.alpha == 0;
    eff_cf = results.eff(cond_cf);
    eff = results.eff(cond);
    DELTAeff = eff_cf-eff;
    NTU = results.NTU(cond);
    data = sortrows([NTU DELTAeff],1);
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of NTU on Effectiveness for Various C$_{r}$ (alpha=0.2)','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(C_r_index(:)),'Location','northwest')
title(legend,'C_r')
set(gca, 'XScale', 'log')
axis([1 100 -0.02 0.1])
saveas(figure(3),['figures_eff','/','NTU for various C_r','.png'],'png');  

%fixed C_r, various alpha
figure;
C_r = 0.2;
for i=1:length(alpha_index)
    cond = results.alpha == alpha_index(i) & results.C_r == C_r;
    cond_cf = results.alpha == 0 & results.C_r == C_r;
    eff_cf = results.eff(cond_cf);
    eff = results.eff(cond);
    DELTAeff = eff_cf-eff;
    NTU = results.NTU(cond);
    data = sortrows([NTU DELTAeff],1);
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of NTU on Effectiveness for Various $\alpha$ ($C_r$=0.2)','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)),'Location','northwest')
title(legend,'alpha')
set(gca, 'XScale', 'log')
axis([1 100 -0.02 0.1])
saveas(figure(4),['figures_eff','/','NTU for various alpha C_r 0.2','.png'],'png');  

%% plotting the effect of C_r

%fixed alpha, various NTU
figure;
alpha = 0.2;
for i=1:length(NTU_index)
    cond = results.NTU == NTU_index(i) & results.alpha == alpha;
    cond_cf = results.NTU == NTU_index(i) & results.alpha == 0;
    eff_cf = results.eff(cond_cf);
    eff = results.eff(cond);
    DELTAeff = eff_cf-eff;
    C_r = results.C_r(cond);
    plot(C_r,DELTAeff);
    hold on
end
title('Effect of Capacitance Ratio on Effectiveness for Various NTU (alpha=0.2)','interpreter','latex')
xlabel('Capacitance Ratio C$_{r}$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)),'Location','northwest')
title(legend,'NTU')
axis([0 1 -0.02 0.1])
saveas(figure(5),['figures_eff','/','C_r for various NTU','.png'],'png');  

%fixed NTU, various alpha
figure;
NTU = 6.3096;
for i=1:length(alpha_index)
    cond = results.alpha == alpha_index(i) & results.NTU == NTU;
    cond_cf = results.alpha == 0 & results.NTU == NTU;
    eff_cf = results.eff(cond_cf);
    eff = results.eff(cond);
    DELTAeff = eff_cf-eff;
    C_r = results.C_r(cond);
    plot(C_r,DELTAeff);
    hold on
end
title('Effect of Capacitance Ratio on Effectiveness for Various $\alpha$ (NTU=6.3096)','interpreter','latex')
xlabel('Capacitance Ratio C$_{r}$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)),'Location','northwest')
title(legend,'alpha')
axis([0 1 -0.02 0.1])
saveas(figure(6),['figures_eff','/','C_r for various alpha','.png'],'png');  

%% compare MATLAB model to theoretical calculations

%counter-flow case
NTU = logspace(0,2,11);
C_r = [0.2 0.5 0.8 1];
figure;
for j=1:length(C_r)
    eff_th = [];
    for i=1:length(NTU)
        if C_r(j) == 1
            eff = NTU(i)/(1+NTU(i));
            eff_th = [eff_th eff];
        else
        eff = (1-exp(-NTU(i)*(1-C_r(j))))/(1-C_r(j)*exp(-NTU(i)*(1-C_r(j))));
        eff_th = [eff_th eff];
        end
    end
    plot(NTU,eff_th,'--','DisplayName',string(C_r(j)))
    hold on
    cond = results.alpha == 0 & results.C_r == C_r(j);
    eff_m = results.eff(cond);
    NTU_m = results.NTU(cond);
    data = sortrows([NTU_m eff_m],1);
    plot(data(:,1),data(:,2),'DisplayName',string(C_r(j)));
end
title('Comparison of MATLAB and Theoretical Effectiveness for Counter-flow','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness $\left[ - \right]$','interpreter','latex')
grid on
text(60,0.7,'- MATLAB')
text(60,0.66,'-- Theory')
axis([1 100 0.5 1.1])
legend('Location','southeast')
title(legend,'C_r')
saveas(figure(7),['figures_eff','/','theory vs matlab counter','.png'],'png');  


%cross-flow case
NTU = logspace(0,2,11);
C_r = [0.2 0.5 0.8 1];
figure;
for j=1:length(C_r)
    eff_th = [];
    for i=1:length(NTU)
        eff = 1-exp(NTU(i)^(0.22)/C_r(j)*(exp(-C_r(j)*NTU(i)^(0.78))-1));
        eff_th = [eff_th eff];
    end
    plot(NTU,eff_th,'--','DisplayName',string(C_r(j)))
    hold on
    cond = results.alpha == 1 & results.C_r == C_r(j);
    eff_m = results.eff(cond);
    NTU_m = results.NTU(cond);
    data = sortrows([NTU_m eff_m],1);
    plot(data(:,1),data(:,2),'DisplayName',string(C_r(j)));
%     figure;
%     hold on
%     plot(eff_th,eff_m,'r');
%     plot(linspace(0,1),linspace(0,1),'k--');
%     plot(linspace(0,0.95*1),linspace(0,1),'k:');
%     plot(linspace(0,1.05*1),linspace(0,1),'k:');
end
title('Comparison of MATLAB and Theoretical Effectiveness for Cross-flow','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness $\left[ - \right]$','interpreter','latex')
grid on
text(60,0.7,'- MATLAB')
text(60,0.66,'-- Theory')
axis([1 100 0.5 1.1])
legend('Location','southeast')
title(legend,'C_r')
saveas(figure(8),['figures_eff','/','theory vs matlab cross','.png'],'png');  

