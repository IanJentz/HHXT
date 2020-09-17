function v = videoFiberPlot(obj)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here

v = VideoWriter([obj.wdir,'\videoFiberPlot.avi']);
open(v);

hfig = obj.initFiberPlot;
frame = getframe(hfig);
writeVideo(v,frame);

times = obj.Fiber.elapsedSeconds'-obj.sec_diff;

for k = 1:length(times)
    hfig = obj.updateFiberPlot(hfig,times(k));
    frame = getframe(hfig);
    writeVideo(v,frame);
end

close(v);

end

