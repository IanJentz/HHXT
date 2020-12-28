cind = model.Runtime.ConvergenceData(:,2) == 1;
cdata = model.Runtime.ConvergenceData(cind,:);
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

% subplot(3,1,1)
plot(steps,cdata(:,3:9)); set(gca,'YScale','log')
% subplot(3,1,2)
% plot(idplot,dcmean(idplot));
% subplot(3,1,3)
% plot(idplot,ddcmean(idplot));

xlabel('number of full steps','Interpreter','latex')
ylabel('convergence value','Interpreter','latex')
lgd = legend('Residual','Abs. $\Delta$: Body $T$','Abs. $\Delta$: Fluid $P$','Abs. $\Delta$: Fluid $T$',...
    'Rel. $\Delta$: Body $T$','Rel. $\Delta$: Fluid $P$','Rel. $\Delta$: Fluid $T$');
 lgd.Interpreter = 'latex';
 lgd.NumColumns = 1;