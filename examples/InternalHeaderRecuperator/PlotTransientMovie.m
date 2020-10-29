
Plot = hhxt.PlotHHXT(filename);

types = {'Temperatures','DarcyVelocityMagnitude','ReynoldsMagnitude','FluidHeating'};
type = types{1};

times = results.SolutionTimes;
t_frames = 1:length(times);

switch type
    case 'ReynoldsMagnitude'
        Cmin = min(mean(results.Reynolds(results.Mesh.findNodes('Region','Face',6),1,1)),mean(results.Reynolds(results.Mesh.findNodes('Region','Face',6),2,1)));
        Cmax = max(mean(results.Reynolds(results.Mesh.findNodes('Region','Face',6),1,end)),mean(results.Reynolds(results.Mesh.findNodes('Region','Face',6),2,end)));
    case 'FluidHeating'
        Cmin = min(min(min(results.StreamVolHeating)));
        Cmax = max(max(max(results.StreamVolHeating)));
        
end

movie_name = [type,'_transient'];
v = VideoWriter(movie_name);
open(v)

for t = t_frames
    
   switch type 
       case 'Temperatures'
           Plot.updatePlot('Temperatures','Mesh','on','Time',t,'TempLimits',[0,T_H_in]);
       case 'DarcyVelocityMagnitude'
           Plot.updatePlot('DarcyVelocityMagnitude','Mesh','on','Time',t);
       case 'ReynoldsMagnitude'
           Plot.updatePlot('ReynoldsMagnitude','Mesh','on','Time',t,'CLimits',[Cmin,Cmax]);
       case 'FluidHeating'
           Plot.updatePlot('FluidHeating','Mesh','on','Time',t,'CLimits',[Cmin,Cmax]);

   end
   
   delete(findall(gcf,'type','annotation'))
   switch type
       
       case 'Temperatures'
        annotation('textbox',[0.6,0.28,0.1,0.1],'String',['time = ',num2str(times(t)),' sec']);
       case {'DarcyVelocityMagnitude','ReynoldsMagnitude','FluidHeating'}
       annotation('textbox',[0.45,0.45,0.1,0.1],'String',['time = ',num2str(times(t)),' sec']);    
   end
   
   frame = getframe(Plot.fighandle);
   writeVideo(v,frame); 
    
    
end



close(v);



