%% Side headers data analysis
%load the results into the workspace
results = readtable('SideHX_short.txt');
close all
%since there are three DOF, results X will be plotted for one variable A
%fixed at a value close to the original case, and the second variable B
%varying for different values of the third variable C
%e.g. plot ineff (X) as a function of alpha (A) for 3 values of C_r (C)
%with NTU (B) fixed at 5
%values for which variable C will plotted
C_r_index = [0.2 0.5 0.7 0.9];
NTU_index = [1 5.2481 18.1970 41.687 63.096];
alpha_index = [0.1 0.2 0.3 0.5];

%% compare MATLAB model to theoretical calculations
color = char('r','b','g','m');
%cross-flow case
NTU = logspace(0,1.8,11);
C_r = [0.2 0.5 0.8 1];
figure;
for j=1:length(C_r)
    eff_th = [];
    for i=1:length(NTU)
        eff = 1-exp(NTU(i)^(0.22)/C_r(j)*(exp(-C_r(j)*NTU(i)^(0.78))-1));
        eff_th = [eff_th eff];
    end
    plot(NTU,eff_th,strcat('-',color(j)),'DisplayName',string(C_r(j)))
    hold on
    cond = results.alpha == 1 & results.C_r == C_r(j);
    eff_m = results.eff(cond);
    NTU_m = results.NTU(cond);
    data = sortrows([NTU_m eff_m],1);
    plot(data(:,1),data(:,2),strcat('*',color(j)),'DisplayName',string(C_r(j)));
end
title('Comparison of MATLAB and Theoretical Effectiveness for Cross-flow','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness $\left[ - \right]$','interpreter','latex')
grid on
text(10,0.7,'* MATLAB')
text(10,0.66,'- Theory')
axis([1 70 0.5 1.1])
legend('Location','southeast')
title(legend,'C_r')
set(gca, 'XScale', 'log')
saveas(figure(1),['figures_delta','/','theory vs matlab cross','.png'],'png');  

%counter-flow case
NTU = logspace(0,1.8,11);
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
    plot(NTU,eff_th,strcat('-',color(j)),'DisplayName',string(C_r(j)))
    hold on
    cond = results.alpha == 0 & results.C_r == C_r(j);
    eff_m = results.eff(cond);
    NTU_m = results.NTU(cond);
    data = sortrows([NTU_m eff_m],1);
    plot(data(:,1),data(:,2),strcat('*',color(j)),'DisplayName',string(C_r(j)));
end
title('Comparison of MATLAB and Theoretical Effectiveness for Counter-flow','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness $\left[ - \right]$','interpreter','latex')
grid on
text(10,0.7,'* MATLAB')
text(10,0.66,'- Theory')
axis([1 70 0.5 1.1])
legend('Location','southeast')
title(legend,'C_r')
set(gca, 'XScale', 'log')
saveas(figure(2),['figures_delta','/','theory vs matlab counter','.png'],'png');  

%% plotting the effect of Area ratio

%fixed NTU, various C_r
figure;
NTU = 5.2481;
for i=1:length(C_r_index)
    cond = results.C_r == C_r_index(i) & results.NTU == NTU;
    [row,col] = find(cond>0);
    eff = results.eff(row);
    if C_r_index(i) == 1
        eff_th = NTU/(1+NTU);
    else
        eff_th = (1-exp(-NTU*(1-C_r_index(i))))/(1-C_r_index(i)*exp(-NTU*(1-C_r_index(i))));
    end
    DELTAeff = eff_th-eff;
    alpha = results.alpha(row);
    plot(alpha,DELTAeff);
    hold on
end
title('Effect of Area Ratio on Effectiveness for Various C$_{r}$ (NTU=5.2481)','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(C_r_index(:)),'Location','northwest')
title(legend,'C_r')
axis([0 1 -0.01 0.1])
saveas(figure(3),['figures_delta','/','alpha for various C_r','.png'],'png');  
%%
%fixed C_r, various NTU
figure;
C_r = 1;
for i=1:length(NTU_index)
    data = [];
    cond = results.NTU == NTU_index(i) & results.C_r == C_r;
    [row,col] = find(cond>0);
    eff = results.eff(row);
    if C_r == 1
        eff_th = NTU_index(i)/(1+NTU_index(i));
    else
        eff_th = (1-exp(-NTU_index(i)*(1-C_r)))/(1-C_r*exp(-NTU_index(i)*(1-C_r)));
    end
    DELTAeff = eff_th-eff;
    alpha = results.alpha(row);
    data = sortrows([alpha DELTAeff],1);
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of Area ratio on Effectiveness for Various NTU ($C_r$=1)','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)),'Location','northwest')
title(legend,'NTU')
axis([0 1 -0.01 0.1])
saveas(figure(4),['figures_delta','/','alpha for various NTU','.png'],'png');  

%% plotting the effect of NTU

%fixed alpha, various C_r
figure;
alpha = 0.2;
for i=1:length(C_r_index)
    data = [];
    cond = results.C_r == C_r_index(i) & results.alpha == alpha;
    [row,col] = find(cond>0);
    eff = results.eff(row);
    NTU = logspace(0,1.8,11);
    data_th =[];
    for j=1:length(NTU)
        if C_r_index(i) == 1
            eff_th = NTU(j)/(1+NTU(j));
        else
            eff_th = (1-exp(-NTU(j)*(1-C_r_index(i))))/(1-C_r_index(i)*exp(-NTU(j)*(1-C_r_index(i))));
        end
        data_th = [data_th ;eff_th];
    end
    NTU = results.NTU(row);
    data = sortrows([NTU eff],1);
    data = [data(:,1) data_th-data(:,2)];
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
axis([1 70 -0.01 0.03])
saveas(figure(5),['figures_delta','/','NTU for various C_r alpha 0.2','.png'],'png');  
%%
%fixed C_r, various alpha
alpha_index = [0.2 0.3 0.4 0.5];
figure;
C_r = 0.9;
for i=1:length(alpha_index)
    data = [];
    cond = results.alpha == alpha_index(i) & results.C_r == C_r;
    [row,col] = find(cond>0);
    eff = results.eff(row);
    NTU = results.NTU(row);
    NTU_th = logspace(0,1.8,11);
    data_th =[];
    for j=1:length(NTU_th)
        if C_r == 1
            eff_th = NTU_th(j)/(1+NTU_th(j));
        else
            eff_th = (1-exp(-NTU_th(j)*(1-C_r)))/(1-C_r*exp(-NTU_th(j)*(1-C_r)));
        end
        data_th = [data_th ;eff_th];
    end
    data = sortrows([NTU eff],1);
    DELTA = data_th-data(:,2);
    data = [data(:,1) DELTA];
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of NTU on Effectiveness for Various $\alpha$ ($C_r$=0.9)','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)),'Location','northwest')
title(legend,'alpha')
set(gca, 'XScale', 'log')
axis([1 70 -0.01 0.05])
saveas(figure(6),['figures_delta','/','NTU for various alpha C_r 0.9','.png'],'png');  

%% plotting the effect of C_r

%fixed alpha, various NTU
figure;
alpha = 0.2;
for i=1:length(NTU_index)
    cond = results.NTU == NTU_index(i) & results.alpha == alpha;
    [row,col] = find(cond>0);
    eff = results.eff(row);
    C_r = results.C_r(row);
    data_th =[];
    for j=1:length(C_r)
        if C_r(j) == 1
            eff_th = NTU_index(i)/(1+NTU_index(i));
        else
            eff_th = (1-exp(-NTU_index(i)*(1-C_r(j))))/(1-C_r(j)*exp(-NTU_index(i)*(1-C_r(j))));
        end
        data_th = [data_th; eff_th];
    end
    data = sortrows([C_r eff],1);
    data = [data(:,1) data_th-data(:,2)];
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of Capacitance Ratio on Effectiveness for Various NTU (alpha=0.2)','interpreter','latex')
xlabel('Capacitance Ratio C$_{r}$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)),'Location','northwest')
title(legend,'NTU')
axis([0 1 -0.01 0.02])
saveas(figure(7),['figures_delta','/','C_r for various NTU','.png'],'png');  
%%
%fixed NTU, various alpha
alpha_index = [0.1 0.2 0.3 0.5];
figure;
NTU = 5.2481;
for i=1:length(alpha_index)
    cond = results.alpha == alpha_index(i) & results.NTU == NTU;
    [row,col] = find(cond>0);
    eff = results.eff(row);
    C_r = results.C_r(row);
    data_th =[];
    for j=1:length(C_r)
        if C_r(j) == 1
            eff_th = NTU/(1+NTU);
        else
            eff_th = (1-exp(-NTU*(1-C_r(j))))/(1-C_r(j)*exp(-NTU*(1-C_r(j))));
        end
        data_th = [data_th; eff_th];
    end
    data = sortrows([C_r eff],1);
    data = [data(:,1) data_th-data(:,2)];
    plot(data(:,1),data(:,2));
    hold on
end
title('Effect of Capacitance Ratio on Effectiveness for Various $\alpha$ (NTU=5.2481)','interpreter','latex')
xlabel('Capacitance Ratio C$_{r}$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)),'Location','northwest')
title(legend,'alpha')
axis([0 1 -0.01 0.04])
saveas(figure(8),['figures_delta','/','C_r for various alpha','.png'],'png');  

