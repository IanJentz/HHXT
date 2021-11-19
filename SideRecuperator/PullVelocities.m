function vel_interps = PullVelocities(model,results,x_corn,y_corn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin <= 2 
x_corn = 0.1836;
y_corn = 0.04575;
end

Nx = 21;
Ny = 51;
vel_interps = cell(1,2);
vel_interps{3} = zeros(Ny,Nx);
vel_interps{4} = zeros(Ny,Nx);

xp = linspace(-x_corn,x_corn,Nx);
yq = linspace(-y_corn,y_corn,Ny);

[vel_interps{1},vel_interps{2}] = meshgrid(xp,yq);


for i = 1:Nx
    
    xq = xp(i)*ones(1,length(yq));
    vel = evaluate(results.InterpolantChannelVelocity,xq,yq);
    
    vel_interps{3}(:,i) = vel(1:Ny);
    vel_interps{4}(:,i) = vel((Ny+1):end);
    
end


end

