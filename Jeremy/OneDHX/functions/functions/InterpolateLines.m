function T_centerline = InterpolateLines(results,L)
% InterpolateLines - Interpolates the HHXT model result at the centerline 
% y = 0.  Outputs results as a table T_centerline

yq = ones(1,50)*0; % the centerline is located at y = 0
xq = L*linspace(-0.5,0.5,50);

% core temp (ind = 1)
ind = 1;
T_core = interpolateSolution(results,xq,yq,ind);

% cold temp (ind = 3)
ind = 3;
T_cold = interpolateSolution(results,xq,yq,ind);

% core temp (ind = 5)
ind =5;
T_hot = interpolateSolution(results,xq,yq,ind);

% cold Press (ind = 2)
ind = 2;
P_cold = interpolateSolution(results,xq,yq,ind);

% hot press (ind = 4)
ind = 4;
P_hot = interpolateSolution(results,xq,yq,ind);

% velocities
pq = [xq;yq];
XVels = evaluate(results.InterpolantXChannelVelocity,pq);
XVels = reshape(XVels,50,2);
XVel_C = XVels(:,1);
XVel_H = XVels(:,2);

T_centerline = [T_hot,T_cold,T_core];
T_centerline = [T_centerline,abs(XVel_C),XVel_H];
T_centerline = [T_centerline,P_cold,P_hot];
T_centerline = [linspace(0,L,50)',T_centerline];

T_centerline = table(T_centerline(:,1),T_centerline(:,2),T_centerline(:,3),T_centerline(:,4),T_centerline(:,5),T_centerline(:,6),T_centerline(:,7),T_centerline(:,8),'VariableNames',{'x','T_H','T_C','T_s','v_C','v_H','P_C','P_H'});

end
