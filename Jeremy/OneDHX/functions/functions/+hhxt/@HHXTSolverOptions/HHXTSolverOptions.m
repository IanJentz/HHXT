classdef (Sealed) HHXTSolverOptions  < handle & matlab.mixin.internal.CompactDisplay 
% HHXTSolverOptions Algorithm options for the HHXT solvers
%
%
%    Instances of this class cannot be created directly, an instance
%    of this class is created by the PDEModel class to provide controls
%    for the solver options.
%
% HHXTSolverOptions properties:
%    AbsoluteTolerance - Absolute error tolerance
%    RelativeTolerance - Relative error tolerance 
%    ResidualTolerance - Acceptable residual tolerance
%    MaxIterations     - Maximum solver iterations
%    MinStep           - Minimum damping of search direction
%    ResidualNorm      - Residual norm
%    ReportStatistics  - Toggle for the statistics report from the solver
%
% See also pde.PDEModel, pde.PDEModel/SolverOptions

% Copyright 2015-2016 The MathWorks, Inc.



properties
% AbsoluteTolerance - Absolute error tolerance
% AbsoluteTolerance is a threshold below which the value of the solution 
% component is unimportant. The absolute error tolerance determines the 
% accuracy when the solution approaches zero.
     AbsoluteTolerance;

% RelativeTolerance - Relative error tolerance 
% This tolerance is a measure of the error relative to the size of each 
% solution component. Roughly, it controls the number of correct digits in 
% all solution components, except those smaller than thresholds imposed by
% the AbsoluteTolerance.  The default, 1e-3, corresponds to 0.1% accuracy.
     RelativeTolerance;

     
% ResidualTolerance - Acceptable residual tolerance
% Residual size at termination of nonlinear solver, specified as a positive 
% scalar. The nonlinear solver iterates until the residual size is less than 
% ResidualTolerance
     ResidualTolerance;

% MaxIterations - Maximum solver iterations
% Maximum number of Gauss-Newton iterations allowed for the nonlinear solver, 
% specified as a positive integer.
     MaxIterations;

% MinStep - Minimum damping of search direction
% Minimum damping of search direction for the nonlinear solver, specified 
% as a positive scalar.
     MinStep;

% ResidualNorm - Residual norm
% Residual norm, specified as the p value for Lp norm, or as the string 'energy'. 
% For the Lp norm, p can be any positive real value, Inf, or -Inf. 
% The p norm of a vector v is sum(abs(v)^p)^(1/p). See norm
     ResidualNorm;
     
     
     
% ReportStatistics - Toggle for the statistics report from the solver
% Print the statistics/convergence information, specified as 'off' or 'on'.
     ReportStatistics 
     
     PlotStatistics
     
     PermeabilityMode
     
     InitializeMode
     
     ResistanceMode
     
     AnisoConductivity
     
     
     
     
end

methods(Hidden=true, Access={?pde.EquationModel})
    function obj=PDESolverOptions()     
      obj.AbsoluteTolerance = 1e-6;
      obj.RelativeTolerance = 1e-3;
      obj.ResidualTolerance = 1e-4; 
      obj.MaxIterations = 25;
      obj.MinStep = 1/2^16;   
      obj.ResidualNorm =  Inf;    
      obj.ReportStatistics  = 'off';         
    end
    function delete(obj)%#ok
    end       
end

methods    
     function set.AbsoluteTolerance(self, atol)     
        self.ValidateAbsoluteTolerance(atol);     
        self.AbsoluteTolerance = atol;       
    end
    
    function set.RelativeTolerance(self, rtol)      
        self.ValidateRelativeTolerance(rtol);      
        self.RelativeTolerance = rtol;        
    end        
    
    function set.ResidualTolerance(self, tol)    
      self.ValidateResidualTolerance(tol);     
      self.ResidualTolerance = tol;     
    end
    
    function set.MaxIterations(self, maxiter) 
      self.ValidateMaxIterations(maxiter);
      self.MaxIterations = maxiter;   
    end
    
    function set.MinStep(self, minstep)  
      self.ValidateMinStep(minstep);      
      self.MinStep = minstep;   
    end
    
    function set.ResidualNorm(self, normval)    
      self.ValidateResidualNorm(normval);  
      if ischar(normval)
           self.ResidualNorm = 'energy';                
      else
        self.ResidualNorm = normval;     
      end
    end
    
    function set.ReportStatistics(self, reportstats)  
      self.ValidateReportStatistics(reportstats);
       nc = numel(reportstats);        
       if strncmpi(reportstats,'on',nc)                  
          self.ReportStatistics = 'on';  
       else
           self.ReportStatistics = 'off';  
       end
    end

end

methods(Static, Access = private)
    function ValidateAbsoluteTolerance(atol)   
      if ~pde.PDESolverOptions.ToleranceOK(atol)
        error(message('pde:pdeModel:invalidAbsoluteTolerance'));   
      end                    
    end
    function ValidateRelativeTolerance(rtol)   
      if ~pde.PDESolverOptions.ToleranceOK(rtol)
        error(message('pde:pdeModel:invalidRelativeTolerance'));   
      end           
    end  
    
    function ValidateResidualTolerance(tol)
      if ~pde.PDESolverOptions.ToleranceOK(tol)
        error(message('pde:pdeModel:invalidResidualTolerance'));   
      end         
    end  
                
    function ValidateMaxIterations(maxiter)    
      if ~(isnumeric(maxiter) && isreal(maxiter) && isscalar(maxiter) && isfinite(maxiter)) ||  maxiter < 0 || issparse(maxiter) || mod(maxiter,1)~=0 
       error(message('pde:pdeModel:invalidMaxIterations'));   
      end           
    end  
   
    function ValidateMinStep(minstep)
      if ~pde.PDESolverOptions.ToleranceOK(minstep)
        error(message('pde:pdeModel:invalidMinStep'));   
      end        
    end  
    
    function ValidateResidualNorm(normval)
       if isnumeric(normval)
            if ~(isnumeric(normval) && isreal(normval) && isscalar(normval))  || issparse(normval) || isnan(normval)
                error(message('pde:pdeModel:invalidResidualNorm'));              
            elseif isfinite(normval)
                if normval < 0
                    error(message('pde:pdeModel:invalidResidualNorm'));  
                end
            end
       elseif ischar(normval)
            nc = numel(normval);        
            if ~strncmpi(normval,'energy',nc) 
                error(message('pde:pdeModel:invalidResidualNorm'))
            end    
       else
           error(message('pde:pdeModel:invalidResidualNorm'))
       end            
    end  
    
    function ValidateReportStatistics(reportstats)     
       if isnumeric(reportstats)           
             error(message('pde:pdeModel:invalidReportStatistics'));              
       elseif ischar(reportstats)
            nc = numel(reportstats); 
            if nc == 1 && strncmpi(reportstats,'o',nc)
                error(message('pde:pdeModel:ambiguousReportStatistics'))
            end
            if ~(strncmpi(reportstats,'on',nc) || strncmpi(reportstats,'off',nc))
                error(message('pde:pdeModel:invalidReportStatistics'))
            end       
       else
           error(message('pde:pdeModel:invalidReportStatistics'))
       end     
    end  
    function ok = ToleranceOK(tol)
      ok = true;
      if ~(isnumeric(tol) && isreal(tol) && isscalar(tol) && isfinite(tol)) ||  tol < 0 || issparse(tol) 
        ok = false;
      end      
    end
    
end    
  
end

