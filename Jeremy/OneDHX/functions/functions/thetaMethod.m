function varargout = thetaMethod(ode,tspan,y0,options,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
varargin = varargin{:};
if nargin == 5
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'theta'
                theta = varargin{i+1};
                varargin(i+1) = [];
                varargin(i) = [];
        end
    end
else
    theta = 0.5; % Crank Nicholson method.
end

solver_name = 'CrankNicholson';
solver_name = 'ode15i';

%update function
nlfcn = options.NonLinearUpdate;

if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error(message('MATLAB:ode15s:NotEnoughInputs'));
      end  
    end
  end
end

% Stats
nsteps   = 0;
nfailed  = 0;
nfevals  = 0; 
npds     = 0;
ndecomps = 0;
nsolves  = 0;

% Output
FcnHandlesUsed  = isa(ode,'function_handle');
output_sol = (FcnHandlesUsed && (nargout==1));      % sol = odeXX(...)
output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)

sol = []; kvec = []; dif3d = []; 
if output_sol
  sol.solver = solver_name;
  sol.extdata.odefun = ode;
  sol.extdata.options = options;                       
  sol.extdata.varargin = varargin;  
end  

% % Handle solver arguments
% t0 = tspan(1); tfinal = tspan(end);
% tdir = sign(tspan(end)-tspan(1));
% hmax = tspan(2)-tspan(1);
% odeFcn = ode;
% neq = size(y0,1);
% ntspan = size(tspan,2);
%
% f0 = feval(odeFcn,t0,y0);
% nfevals = nfevals + 1;

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
 options, threshold, rtol, normcontrol, normy, hmax, htry, htspan] = ...
    odearguments(FcnHandlesUsed, 'ode15s', ode, tspan, y0, options, varargin);
nfevals = nfevals + 1;
one2neq = (1:neq);

% Handle the mass matrix
[Mtype, Mt, Mfun, Margs, dMoptions] = odemass(FcnHandlesUsed,odeFcn,t0,y0,...
                                              options,varargin);


                                          
refine = max(1,odeget(options,'Refine',1,'fast'));
if ntspan > 2
  outputAt = 'RequestedPoints';         % output only at tspan points
elseif refine <= 1
  outputAt = 'SolverSteps';             % computed points, no refinement
else
  outputAt = 'RefinedSteps';            % computed points, with refinement
  S = (1:refine-1) / refine;
end
printstats = strcmp(odeget(options,'Stats','off','fast'),'on');

% Handle the Jacobian
[Jconstant,Jac,Jargs,Joptions] = ...
    odejacobian(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
Janalytic = isempty(Joptions);

t = t0; y = y0; 
f = f0;

updated = feval(nlfcn,t,y); % update the permeability


Mcurrent = true;
Mtnew = Mt;

% Adjust the warnings.
warnoffId = { 'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'}; 
for i = 1:length(warnoffId)    
  warnstat(i) = warning('query',warnoffId{i});
  warnoff(i) = warnstat(i);
  warnoff(i).state = 'off';
end

if Jconstant 
  if isempty(Jac) % use odenumjac
    [Jac,Joptions.fac,nF] = odenumjac(odeFcn, {t0,y0,odeArgs{:}}, f0, Joptions);    
    nfevals = nfevals + nF;
    npds = npds + 1;
  elseif ~isa(Jac,'numeric')  % not been set via 'options'  
    Jac = feval(Jac,t0,y0,Jargs{:}); % replace by its value
    npds = npds + 1;
  end
end

    if Jconstant
        dfdy = Jac;
    elseif Janalytic
        dfdy = feval(Jac,t,y,Jargs{:});     
        npds = npds + 1;                            
    else   % Joptions not empty
        [dfdy,Joptions.fac,nF] = odenumjac(odeFcn, {t,y,odeArgs{:}}, f0, Joptions);  
        nfevals = nfevals + nF;    
        npds = npds + 1;                            
    end   
Jcurrent = true;

maxk = 1; % maximum order of time differentiation is 1st order


absh = options.AbsTol/max(abs(f0));  %initial guess for step size


% Adjust the warnings.
warnoffId = { 'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'}; 
for i = 1:length(warnoffId)    
  warnstat(i) = warning('query',warnoffId{i});
  warnoff(i) = warnstat(i);
  warnoff(i).state = 'off';
end

% Initialize.
k = 1;                                  % start at order 1 with BDF1
K = 1;                                  % K = 1:k
klast = k;
abshlast = absh;

invGa(k) = 1;

hinvGk = absh*tdir * invGa(k);

Miter = Mt - hinvGk * dfdy;
Miternew = Miter;

if issparse(Miter)
    [L,U,P,Q,R] = lu(Miter);
    Lnew = L; Unew = U; Pnew = P; Qnew = Q; Rnew = Q;
else
    [L,U,p] = lu(Miter,'vector');
    Lnew = L; Unew = U; pnew = p;
end
ndecomps = ndecomps + 1;


dif = zeros(neq,maxk+2);
dif(:,1) = absh * f0;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
  if output_sol
    chunk = min(max(100,50*refine), refine+floor((2^11)/neq));      
    tout = zeros(1,chunk);
    yout = zeros(neq,chunk);
    kvec = zeros(1,chunk);
    dif3d = zeros(neq,maxk+2,chunk);
  else      
    if ntspan > 2                         % output only at tspan points
      tout = zeros(1,ntspan);
      yout = zeros(neq,ntspan);
    else                                  % alloc in chunks
      chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
      tout = zeros(1,chunk);
      yout = zeros(neq,chunk);
    end
  end  
  nout = 1;
  tout(nout) = t;
  yout(:,nout) = y;  
end




%% THE MAIN LOOP
absh = 0.1;

iter = 1;
done = false;
at_hmin = false;
while ~done
  
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));

    h = tdir*absh;
    
    % Stretch the step if within 10% of tfinal-t.
    if 1.1*absh >= abs(tfinal - t)
        h = tfinal - t;
        absh = abs(h);
        done = true;
    end

    recalcMiter = false;
    % LOOP FOR ADVANCING ONE STEP
    while true
        % Predict a solution at t+h.
        tnew = t + h;
        if done
            tnew = tfinal;   % Hit end point exactly.
        end
        h = tnew - t;      % Purify h.
        hinvGk = h * invGa(k);
        
        f = feval(odeFcn,t,y);
        nfevals = nfevals + 1;
        
       

        Miter = Mt - hinvGk * dfdy;

        [lastmsg,lastid] = lastwarn('');
        warning(warnoff);
        if issparse(Miter)
          ydt = Q * (U \ (L \ (P * (R \ f))));
        else  
          ydt = U \ (L \ f(p));
        end 
        warning(warnstat);
        
        ypred = y + h*ydt;

        fnew = feval(odeFcn,tnew,ypred);
        nfevals = nfevals + 1;
        
        if Mtype >= 3
            Mtnew = feval(Mfun,tnew,ypred,Margs{:}); % state-dependent
            nfevals = nfevals + 1;
            recalcMiter = true;
        end
        
        if Jconstant
            dfdynew = Jac;
        elseif Janalytic
            dfdynew = feval(Jac,tnew,ypred,Jargs{:});     
            npds = npds + 1; 
            recalcMiter = true;
        else   % Joptions not empty
            [dfdynew,Joptions.fac,nF] = odenumjac(odeFcn, {tnew,ypred,odeArgs{:}}, f, Joptions);  
            nfevals = nfevals + nF;    
            npds = npds + 1;
            recalcMiter = true;
        end  
             
        if recalcMiter
            Miternew = Mtnew - hinvGk * dfdynew;
            if issparse(Miternew)
                [Lnew,Unew,Pnew,Qnew,Rnew] = lu(Miternew);
            else
                [Lnew,Unew,pnew] = lu(Miternew,'vector');
            end
            ndecomps = ndecomps + 1;
            recalcMiter = false;
        end

         
        [lastmsg,lastid] = lastwarn('');
        warning(warnoff);
        if issparse(Miternew)
          ydtnew = Qnew * (Unew \ (Lnew \ (Pnew * (Rnew \ fnew))));
        else  
          ydtnew = Unew \ (Lnew \ fnew(pnew));
        end 
        warning(warnstat);
        
        ydt = theta * ydtnew + (1-theta)* ydt;

        ynew = y + ydt*h;


        hsug = options.AbsTol/max(abs(ydt));


        nsolves = nsolves + iter; 
        
        err = (ynew - ypred)./ynew;
        accl = 0.9;
        
%         if true
%             h = 0.1;
%             if t > 5
%                 h = 1;
%             end
%             break;
%         else

        if true
            
%             range = [ -inf,4 ; 4,6 ; 6,50 ; 50,inf];          
%             hvalu = [ 0.1 ; .05 ; 1.0 ; 5.0];
%             h = hvalu(sum(t <= range,2) == 1);

            range = [ -inf,6 ; 6,29 ; 29,36 ; 36,59 ; 59,inf];          
            hvalu = 0.1*[ 0.1    ; 1.0  ; 0.1   ; 1.0   ; 5.0];
            h = hvalu(sum(t <= range,2) == 1);
            
            break;
        else
            
        if max(err) <= options.AbsTol
            hnew = h*(options.AbsTol/max(err));
            h = (1-accl)*h + accl*hnew + hsug;
            break;
        else
            nfailed = nfailed + 1;
            hnew = h*(options.AbsTol/max(err));
            h = (1-accl)*h + accl*hnew - hsug;
            if h < hsug
                h = hsug;
            end
        end
        
        end
    
    end
    
    nr2 = 1;
    fprintf('%4i%20.4e%12.7f\n',iter,sqrt(nr2),t);
    
    absh = abs(h);
    t = tnew;
    y = ynew;
    dif(:,k) = ydt*h;
    
    
    
    updated = feval(nlfcn,t,y); % update the permeability

    if Mtype >= 3
        Mt = feval(Mfun,t,y,Margs{:}); % state-dependent
        nfevals = nfevals + 1;
        recalcMiter = true;
    end

    if Jconstant
        dfdy = Jac;
    elseif Janalytic
        dfdy = feval(Jac,t,y,Jargs{:});     
        npds = npds + 1; 
        recalcMiter = true;
    else   % Joptions not empty
        [dfdy,Joptions.fac,nF] = odenumjac(odeFcn, {t,y,odeArgs{:}}, f, Joptions);  
        nfevals = nfevals + nF;    
        npds = npds + 1;
        recalcMiter = true;
    end  

    if recalcMiter
        Miter = Mt - hinvGk * dfdy;
        if issparse(Miter)
            [L,U,P,Q,R] = lu(Miter);
        else
            [L,U,p] = lu(Miter,'vector');
        end
        ndecomps = ndecomps + 1;
        recalcMiter = false;
    end
    
    
    nsteps = nsteps + 1;  
    
    if output_sol
        nout = nout + 1;
        if nout > length(tout)
          tout = [tout, zeros(1,chunk)];  % requires chunk >= refine
          yout = [yout, zeros(neq,chunk)];    
          kvec = [kvec, zeros(1,chunk)];
          dif3d = cat(3,dif3d, zeros(neq,maxk+2,chunk));
        end
        tout(nout) = tnew;
        yout(:,nout) = ynew;
        kvec(nout) = k;
        dif3d(:,:,nout) = dif;        
    end   
    
    iter = iter + 1;
%     if ~mod(iter,10)
%         keyboard
%     end

    
end

if output_sol
    sol = struct('solver',solver_name,'extdata',[],'x',tout(1:nout),'y',yout(:,1:nout),'stats',[],'idata',[]);
    sol.extdata = struct('odefun',ode,'options',options,'varargin',{varargin});
    sol.stats = struct('nsteps',nsteps,'nfailed',nfailed,'nfevals',nfevals,'npds',npds,'ndecomps',ndecomps,'nsolves',nsolves);
    sol.idata = struct('kvec',kvec(1:nout),'dif3d',dif3d(:,:,1:nout),'idxNonNegative',[]);
    varargout{1} = sol;
else
    varargout = {};
end


end

