function m=firstOrderODEm(t, u)
%firstOrderODEm - Mass matrix function for first-order ODE system
%Computes the mass matrix of discretized PDE with first-order time
%derivative. 
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2016 The MathWorks, Inc.

global femodel

if(femodel.vd || femodel.vh)
  if(femodel.numConstrainedEqns == femodel.totalNumEqns)
    uFull = u;
  else
    uFull = femodel.B*u + femodel.ud;
  end
end

if(femodel.vd)
        if(femodel.nrp==2)
            femodel.Mass = formGlobalM2D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'d');
        elseif(femodel.nrp==3)
            femodel.Mass = formGlobalM3D(femodel.emptyPDEModel, femodel.p, femodel.t, femodel.coefstruct,uFull,t,'d');
        end    
end

if(femodel.vh)
  [~,~,H,~]=femodel.thePde.assembleBoundary(uFull,t);
  femodel.B=pdenullorth(H);
end
m=femodel.B'*femodel.Mass*femodel.B;

end
