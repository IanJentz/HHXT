function f=firstOrderODEf(t,u)
%firstOrderODEf - Residue function for first-order ODE system
%Computes residue of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel

if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || ...
        femodel.vc || femodel.va || femodel.vf)
    f=-femodel.K*u+femodel.F;
    return
end

if(femodel.numConstrainedEqns == femodel.totalNumEqns)
    uFull = u;
else
    uFull = femodel.B*u + femodel.ud;
end

% if femodel.vc || femodel.va || femodel.vf || femodel.vq || femodel.vg || femodel.vh || femodel.vr
%     if false%t > femodel.thePde.HHXTSolverOptions.curtime
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t;
%     end
% end

% if femodel.vc || femodel.va || femodel.vf || femodel.vq || femodel.vg || femodel.vh || femodel.vr
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t;
% end

if femodel.vq || femodel.vg || femodel.vh || femodel.vr
    [femodel.Q,femodel.G,femodel.H,femodel.R]=femodel.thePde.assembleBoundary(uFull,t);
    
    % deal with outflow boundaries
    femodel.outflow_ind = any(isnan(femodel.Q)')';
    femodel.Q(femodel.outflow_ind,:) = 0;    
end

if femodel.vc || femodel.va || femodel.vf
    
%     if true%t > femodel.thePde.HHXTSolverOptions.curtime
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t;
%     end
    
    if(femodel.nrp==2)
        [K,F] = formGlobalKF2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
        [B,KB, FB, AB] = formGlobalSUPG2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct, femodel.coefstructb,uFull,t);
    elseif(femodel.nrp==3)
        [K,F] = formGlobalKF3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
        [B,KB, FB, AB] = formGlobalSUPG3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct, femodel.coefstructb,uFull,t);
    end
    
    % deal with outflow boundaries
    KB(femodel.outflow_ind,:) = 0;
    B(femodel.outflow_ind,:) = -B(femodel.outflow_ind,:);
    
    % add in petrov corrections
    femodel.K = K + B + KB; % b terms and Petrov corrections
    femodel.A = A + AB; % Petrov corrections
    femodel.F = F + FB; % Petrov corrections          
    
end


if ~(femodel.vh || femodel.vr)
    % neither H nor R are functions of t or u
    K = femodel.K + femodel.A + femodel.Q;
    F = femodel.B'*(femodel.F + femodel.G - K*femodel.ud);
    K = femodel.B'*K*femodel.B;
else
    dt=femodel.tspan*sqrt(eps);
    t1=t+dt;
    dt=t1-t;
    if femodel.vd
        if(femodel.nrp==2)
            femodel.Mass = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'d');
        elseif(femodel.nrp==3)
            femodel.Mass = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'d');
        end
    end
    if femodel.vh
        [N,O]=pdenullorth(femodel.H);
        ud=O*((femodel.H*O)\femodel.R);
        [~,~,H,R]=femodel.thePde.assembleBoundary(uFull,t1);
        [N1,O1]=pdenullorth(H);
        ud1=O1*((H*O1)\R);
        K=N'*((femodel.K+femodel.A+femodel.Q)*N + femodel.Mass*(N1-N)/dt);
        F=N'*((femodel.F+femodel.G)-(femodel.K+femodel.A+femodel.Q)*ud - ...
            femodel.Mass*(ud1-ud)/dt);
    else
        HH=femodel.H*femodel.Or;
        ud=femodel.Or*(HH\femodel.R);
        [~,~,~,R]=femodel.thePde.assembleBoundary(uFull,t1);
        ud1=femodel.Or*(HH\R);
%         K=femodel.Nu'*(femodel.K+femodel.A+femodel.Q)*femodel.Nu;
%         F=femodel.Nu'*((femodel.F+femodel.G)-(femodel.K+femodel.A+femodel.Q)*ud - ...
%             femodel.Mass*(ud1-ud)/dt);
        K=femodel.B'*(femodel.K+femodel.A+femodel.Q)*femodel.B;
        F=femodel.B'*((femodel.F+femodel.G)-(femodel.K+femodel.A+femodel.Q)*ud - ...
            femodel.Mass*(ud1-ud)/dt);
    end
    femodel.ud=ud;
end

f=-K*u+F;

% if femodel.vc || femodel.va || femodel.vf || femodel.vq || femodel.vg || femodel.vh || femodel.vr
%     if t > femodel.thePde.HHXTSolverOptions.curtime
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t;
%     end
% end

end