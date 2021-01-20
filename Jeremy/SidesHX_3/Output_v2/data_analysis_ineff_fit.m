%% Side headers data analysis
%load the results into the workspace
results = readtable('SideHXoutput3_v2.txt');
close all
%since there are three DOF, results X will be plotted for one variable A
%fixed at a value close to the original case, and the second variable B
%varying for different values of the third variable C
%e.g. plot ineff (X) as a function of alpha (A) for 3 values of C_r (C)
%with NTU (B) fixed at 5
%values for which variable C will plotted
NTU_index = [1.5849 6.3096 25.119 39.8110 63.096];
alpha_index = [0.1 0.2 0.3 0.4];
C_r = 0.9;

%% plotting the effect of Area ratio
%fixed C_r, various NTU
figure(1);
figure(2);
h = [];
for i=1:length(NTU_index)
    %plot ineff
    cond = results.NTU == NTU_index(i) & results.C_r == C_r;
    [row,col] = find(cond>0);
    ineff = results.ineff(row);
    alpha = results.alpha(row);
    fit_alpha = fit(alpha,ineff,'poly2');
    figure(1);
    h(i) = plot(alpha,ineff,'*');
    hold on
    plot(fit_alpha)
    %plot eff-loss
    data = [];
    if C_r == 1
        eff_th = NTU_index(i)/(1+NTU_index(i));
    else
        eff_th = (1-exp(-NTU_index(i)*(1-C_r)))/(1-C_r*exp(-NTU_index(i)*(1-C_r)));
    end
    eff = 1-fit_alpha(alpha);
    DELTAeff = eff_th-eff;
    data = sortrows([alpha DELTAeff],1);
    figure(2);
    plot(data(:,1),data(:,2));
    hold on
end
figure(1);
title('Effect of Area ratio on Ineffectivness for Various NTU ($C_r$=0.9)','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
title(legend,'NTU')
legend(h,'1.5136','5.248','18.197','41.687','63.096');
axis([0 1 -0.1 0.5])
saveas(figure(1),['figures_fit','/','ineff as function of alpha C_r 0.9','.png'],'png'); 
figure(2);
title('Effect of Area ratio on Effectiveness for Various NTU ($C_r$=0.9)','interpreter','latex')
xlabel('Area Ratio $\alpha$ $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(NTU_index(:)),'Location','northwest')
title(legend,'NTU')
axis([0 1 -0.01 0.1])
saveas(figure(2),['figures_fit','/','eff loss as function of alpha C_r 0.9','.png'],'png');  

%% plotting the effect of NTU
%fixed C_r, various alpha
figure(3);
figure(4);
h = [];
for i=1:length(alpha_index)
    data = [];
    cond = results.alpha == alpha_index(i) & results.C_r == C_r;
    [row,col] = find(cond>0);
    ineff = results.ineff(row);
    NTU = results.NTU(row);
    data = sortrows([NTU ineff],1);
    fit_NTU = fit(data(:,1),data(:,2),'power2');
    figure(3);
    h(i) = plot(data(:,1),data(:,2),'*');
    hold on
    plot(fit_NTU)
    %plot eff-loss
    data = [];
    data_th =[];
    for j=1:length(NTU)
        if C_r == 1
            eff_th = NTU(j)/(1+NTU(j));
        else
            eff_th = (1-exp(-NTU(j)*(1-C_r)))/(1-C_r*exp(-NTU(j)*(1-C_r)));
        end
        data_th = [data_th ;eff_th];
    end
    eff = 1-fit_NTU(NTU);
    DELTAeff = data_th-eff;
    data = sortrows([NTU DELTAeff],1);
    figure(4);
    plot(data(:,1),data(:,2));
    hold on
end
figure(3);
title('Effect of NTU on Ineffectivness for Various $\alpha$ ($C_r$=0.9)','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Ineffectiveness $\left[ - \right]$','interpreter','latex')
grid on
title(legend,'alpha')
legend(h,'0.1','0.2','0.3','0.4')
set(gca, 'XScale', 'log')
axis([1 70 -0.1 0.5])
saveas(figure(3),['figures_fit','/','ineff as function of NTU C_r 0.9','.png'],'png');  
figure(4);
title('Effect of NTU on Effectiveness for Various $\alpha$ ($C_r$=0.9)','interpreter','latex')
xlabel('NTU $\left[ - \right]$','interpreter','latex')
ylabel('Effectiveness Loss $\left[ - \right]$','interpreter','latex')
grid on
legend(num2str(alpha_index(:)),'Location','northwest')
title(legend,'alpha')
set(gca, 'XScale', 'log')
axis([1 70 -0.01 0.05])
saveas(figure(4),['figures_fit','/','eff loss as function of NTU C_r 0.9','.png'],'png');  

