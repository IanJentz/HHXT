function updated = nonLinearUpdate(t,u)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

updated = false;
global femodel

if(femodel.numConstrainedEqns == femodel.totalNumEqns)
    uFull = u;
else
    uFull = femodel.B*u + femodel.ud;
end

if femodel.vc || femodel.va || femodel.vf || femodel.vq || femodel.vg || femodel.vh || femodel.vr
%     if t > femodel.thePde.HHXTSolverOptions.curtime
    % update the permeability
    tempresult = hhxt.HHXTResults(femodel.thePde,uFull,0);
    femodel.thePde = updatePermeability(femodel.thePde,tempresult);
    femodel.thePde.HHXTSolverOptions.curtime = t;
    updated = true;
%     end
end

end

