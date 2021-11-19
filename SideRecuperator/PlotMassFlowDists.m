function PlotMassFlowDists(filename,dir_plot,breakpoints,c_rho,strm,rid_int)

fNames = {'Cold CO$_2$','Hot CO$_2$'};

sfigfile = false;
mode = 'elemental';

if nargin == 1
    breakpoints = false;
    dir_plot = cd;
end
if nargin == 2
    breakpoints = false;
end
warning off
mkdir([dir_plot,'/figs']);
warning on

load(filename,'results','model','mdot_interps')

[pm,em,tm] = model.Mesh.meshToPet;

H_in = results.Tables{1,4}.edge{1};% 37;
H_out = results.Tables{1,4}.edge{2};%32;
C_in = results.Tables{1,3}.edge{1};%38;
C_out = results.Tables{1,3}.edge{2};%31;

BCs = [C_in,C_out,H_in,H_out];
linestyle = {'ob','xb','or','xr'};

hfig = figure();

hold on

for i = 1:length(BCs)
    
    nodes = model.Mesh.findNodes('region','Edge',BCs(i));
    if i > 2
    dim =2;
    else
    dim = 1;
    end
    x = model.Mesh.Nodes(dim,nodes);
    
    if i <= 2
        temp = results.NodalSolution(nodes,3)';
        press = results.NodalSolution(nodes,2)';
    else
        temp = results.NodalSolution(nodes,5)';
        press = results.NodalSolution(nodes,4)';
    end 
    if isa(c_rho,'function_handle')
        rho = c_rho(1,1,temp,press);
    else
        rho = c_rho*ones(size(temp));
    end
    
    switch mode
        case 'nodal'
        x = (x-min(x))/(max(x)-min(x));
        if i <= 2
        v = results.YChannelVelocity(nodes,1)';
        else
        v = results.YChannelVelocity(nodes,2)';    
        end
        plot(x,-v.*rho,linestyle{i})
        
        case 'elemental'
        elems = find((ismember(tm(1,:),nodes)+ismember(tm(2,:),nodes)+ismember(tm(3,:),nodes)) == 2);
        xe = ( pm(dim,tm(1,elems))+ pm(dim,tm(2,elems))+pm(dim,tm(3,elems)) )/3;     
        x = (xe-min(x))/(max(x)-min(x));
        if i <= 2
        v = results.ElemMassFlux(elems,1)';
        else
        v = results.ElemMassFlux(elems,2)';   
        end 
        temp = 0.5*(temp(1:end-1)+temp(2:end));
        press = 0.5*(press(1:end-1)+press(2:end));
       
        if i == 1; x = 1-x; end
        plot(x,abs(v),linestyle{i})
        
    end
    
    

end

hold off

xlabel('$x^*$','Interpreter','latex')
ylabel('$ j_m \left [ \frac{kg}{m^2s} \right ]$','Interpreter','latex')
hfig.Units = 'inches'; hfig.RendererMode = 'manual'; hfig.Renderer = 'painters';
hfig.OuterPosition = [2,3,3.34,3.75];
hfig.PaperSize = hfig.OuterPosition([end-1,end]);
lgnd = legend('$C_{in}$','$C_{out}$','$H_{in}$','$H_{out}$','Interpreter','latex','Location','northoutside');
lgnd.NumColumns = 2;
% lgnd.Position = [.29,.81,.45,.12];
    
if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

if sfigfile == true; saveas(hfig,[dir_plot,'/figs/BoundaryVelocities.fig'],'fig'); end
saveas(hfig,[dir_plot,'/figs/BoundaryMassFlux.png'],'png');  
saveas(hfig,[dir_plot,'/figs/BoundaryMassFlux.pdf'],'pdf');

close(hfig)

switch mode
    case 'nodal'
    mlims = [ min( min(min(mdot_interps{3})),min(min(mdot_interps{4})) ),...
              max( max(max(mdot_interps{3})),max(max(mdot_interps{4})) ) ];

    hfig = figure();
    subplot(1,2,1)
    contour3(mdot_interps{3},mdot_interps{2},mdot_interps{1},'b');
    view(-90,90);
    xlim(mlims)
    xlabel('$ j_m \left [ \frac{kg}{m^2s} \right ]$','Interpreter','latex')
    ylabel('$y \left [ m \right ]$','Interpreter','latex')
    zlabel('$x \left [ m \right ]$','Interpreter','latex')
    title(fNames{strm(1)},'Interpreter','latex')

    subplot(1,2,2)
    contour3(mdot_interps{4},mdot_interps{2},mdot_interps{1},'r');
    view(-90,90);
    xlim(mlims)
    xlabel('$ j_m \left [ \frac{kg}{m^2s} \right ]$','Interpreter','latex')
    ylabel('$y \left [ m \right ]$','Interpreter','latex')
    zlabel('$x \left [ m \right ]$','Interpreter','latex')
    title(fNames{strm(2)},'Interpreter','latex')
    
    case 'elemental'
    elems = model.Mesh.findElements('region','Face',rid_int);
    x = ( pm(1,tm(1,elems))+pm(1,tm(2,elems))+pm(1,tm(3,elems)) )/3;
    y = ( pm(2,tm(2,elems))+pm(2,tm(2,elems))+pm(2,tm(3,elems)) )/3;
    Nlines = 11;  Npline = 101;
    X = linspace(min(x),max(x),Nlines);
    Y = linspace(min(y),max(y),Npline);
    [X,Y] = meshgrid(X,Y);
    
    vlims = [ min( min(min(results.ElemMassFlux(elems,1))),min(min(results.ElemMassFlux(elems,2))) ),...
              max( max(max(results.ElemMassFlux(elems,1))),max(max(results.ElemMassFlux(elems,2))) ) ];
    
    hfig = figure();
    
    subplot(1,2,1)
    v = results.ElemMassFlux(elems,1)';
    warning off
    F = scatteredInterpolant(x',y',v');
    warning on
    F.Method = 'nearest'; F.ExtrapolationMethod = 'none';
    V = X;
    V = F(X,Y);
    contour3(V,Y,X,'b');
    view(-90,90);
    xlim(vlims)
    xlabel('$ j_m \left [ \frac{kg}{m^2s} \right ]$','Interpreter','latex')
    ylabel('$y \left [ m \right ]$','Interpreter','latex')
    zlabel('$x \left [ m \right ]$','Interpreter','latex')
    title(fNames{strm(1)},'Interpreter','latex')
    
    subplot(1,2,2)
    v = results.ElemMassFlux(elems,2)';
    warning off
    F = scatteredInterpolant(x',y',v');
    warning on
    F.Method = 'nearest'; F.ExtrapolationMethod = 'none';
    V = F(X,Y);
    contour3(V,Y,X,'r');
    view(-90,90);
    xlim(vlims)
    xlabel('$ j_m \left [ \frac{kg}{m^2s} \right ]$','Interpreter','latex')
    ylabel('$y \left [ m \right ]$','Interpreter','latex')
    zlabel('$x \left [ m \right ]$','Interpreter','latex')
    title(fNames{strm(2)},'Interpreter','latex')
        
end

hfig.Units = 'inches'; hfig.RendererMode = 'manual'; hfig.Renderer = 'painters';
hfig.OuterPosition = [2,3,6.85,3.75];
hfig.PaperSize = hfig.OuterPosition([end-1,end]);
    
if breakpoints == true
disp('....press F5 to continue')
keyboard %programatically inserts a breakpoint, comment this out if you don't want breakpoints
end

if sfigfile == true; saveas(hfig,[dir_plot,'/figs/ChannelMassFlow.fig'],'fig'); end
saveas(hfig,[dir_plot,'/figs/ChannelMassFlux.png'],'png');  
saveas(hfig,[dir_plot,'/figs/ChannelMassFlux.pdf'],'pdf');

close(hfig)

end