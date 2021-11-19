function mdot_interps = PullMassFlows(model,results,rho_fun,x_corn,y_corn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 2 
x_corn = 0.1836;
y_corn = 0.04575;
end

Nx = 21;
Ny = 51;
mdot_interps = cell(1,2);
mdot_interps{3} = zeros(Ny,Nx);
mdot_interps{4} = zeros(Ny,Nx);

xp = linspace(-x_corn,x_corn,Nx);
yq = linspace(-y_corn,y_corn,Ny);

[mdot_interps{1},mdot_interps{2}] = meshgrid(xp,yq);


for i = 1:Nx
    
    xq = xp(i)*ones(1,length(yq));
    vel = evaluate(results.InterpolantDarcyFlux,xq,yq);
    T_C = results.interpolateSolution(xq,yq,3);
    P_C = results.interpolateSolution(xq,yq,2);
    T_H = results.interpolateSolution(xq,yq,5);
    P_H = results.interpolateSolution(xq,yq,4);
    if isa(rho_fun,'function_handle')
    rho_C = rho_fun(1,1,T_C,P_C);
    rho_H = rho_fun(1,1,T_H,P_H);
    else
    rho_C = rho_fun*ones(size(T_C));
    rho_H = rho_C;
    end
    
    
    mdot_interps{3}(:,i) = vel(1:Ny).*rho_C;
    mdot_interps{4}(:,i) = vel((Ny+1):end).*rho_H;
    
end




end