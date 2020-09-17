function v = videoMap(obj)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here

v = VideoWriter([obj.wdir,'\videoMap.avi']);
open(v);

hfig = obj.initMap;
frame = getframe(hfig);
writeVideo(v,frame);

times = [obj.Fiber.elapsedSeconds'-obj.sec_diff,obj.seconds];
times = unique(sort(times));

for k = 1:length(times)
    hfig = obj.updateMap(hfig,times(k));
    frame = getframe(hfig);
    writeVideo(v,frame);
end

close(v);

end

