function df=firstOrderODEdf(t,u)
%firstOrderODEdf - Jacobian function for first-order ODE system
%Computes jacobian of discretized PDE with first-order time derivative.
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel

if ~(femodel.vq || femodel.vg || femodel.vh || femodel.vr || femodel.vc || femodel.va || femodel.vf)
  df=-femodel.K;
  return
end

if(femodel.numConstrainedEqns == femodel.totalNumEqns)
  uFull = u;
else
  uFull = femodel.B*u + femodel.ud;
end

% if femodel.vc || femodel.va || femodel.vq || femodel.vh
%     if false%t > femodel.thePde.HHXTSolverOptions.curtime
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t; 
%     end
% end

if(femodel.vq || femodel.vh)
    [Q,~,H,~]=femodel.thePde.assembleBoundary(uFull,t);
    
    % deal with outflow boundaries
    femodel.outflow_ind = any(isnan(Q)')';
    Q(femodel.outflow_ind,:) = 0;
end

if(femodel.vc || femodel.va)
    
%     if t > femodel.thePde.HHXTSolverOptions.curtime
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t;
%     end
    
    
    if(femodel.nrp==2)
        [K, ~] = formGlobalKF2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
        [B,KB, ~, AB] = formGlobalSUPG2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct, femodel.coefstructb,uFull,t);
    elseif(femodel.nrp==3)
        [K, ~] = formGlobalKF3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t);
        A = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'a');
        [B,KB, ~, AB] = formGlobalSUPG3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct, femodel.coefstructb,uFull,t);
    end
    
    % deal with outflow boundaries
    KB(femodel.outflow_ind,:) = 0;
    B(femodel.outflow_ind,:) = -B(femodel.outflow_ind,:);
    
    K = K + B + KB; % b terms and Petrov corrections
    A = A + AB; % Petrov corrections
    
%     % make sure dofs within regions of no definition are handled
%     K = K + femodel.Ni;
    
end

if femodel.vq
  femodel.Q=Q;
end

if femodel.vh
  femodel.B=pdenullorth(H);
end

if femodel.vc
  femodel.K=K;
end

if femodel.va
  femodel.A=A;
end

if ~femodel.vh
  df=-femodel.B'*(femodel.K+femodel.A+femodel.Q)*femodel.B;
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
  [~,~,H,~]=femodel.thePde.assembleBoundary(uFull,t1);
  N1=pdenullorth(H);
  N1 = N1(:,femodel.Bkeep);  % remove 
  df=-femodel.B'*((femodel.K+femodel.A+femodel.Q)*femodel.B + ...
    femodel.Mass*(N1-femodel.B)/dt);
end

% if femodel.vc || femodel.va || femodel.vq || femodel.vh
%     if t > femodel.thePde.HHXTSolverOptions.curtime
%     % update the permeability
%     tempresult = hhxt.HHXTResults(femodel.thePde,uFull);
%     femodel.thePde = updatePermeability(femodel.thePde,tempresult);
%     femodel.thePde.HHXTSolverOptions.curtime = t; 
%     end
% end

end

