classdef (Sealed) CoefficientAssignment  < pde.EquationAssignment & matlab.mixin.internal.CompactDisplay 
% CoefficientAssignment Specify all PDE coefficients over a domain or subdomain
%     The PDE toolbox can solve equations of the form:
%  
%              m*d^2u/dt^2 + d*du/dt - div(c*grad(u)) + a*u = f
%  
%     and the corresponding eigenvalue equation of the form:
%  
%                   - div(c*grad(u)) + a*u = lamda*d*u
%     or
%                   - div(c*grad(u)) + a*u = (lamda^2)*m*u
%  
%     The equation to solve is defined in terms of the coefficients m, d, c, 
%     a, f. This method creates an object representing the coefficients 
%     in a domain or subdomain and appends the object to the EquationCoefficients 
%     property.
%  
%     You can define more than one set of equation coefficients. For example, 
%     for a 2-D geometry that has two subdomains and equation coefficients
%     that differ in each subdomain. However, the same coefficient terms must
%     be present in each equation - they just differ in value.
%     Note: Subdomains are not currently supported in 3-D.
%
%     Instances of this class can only be created by calling the 
%     specifyCoefficients method of the PDEModel class. 
%
% See also pde.PDEModel, pde.PDEModel/specifyCoefficients

% Copyright 2015-2017 The MathWorks, Inc.


properties (SetAccess = private)
% RegionType - Type of geometric region the coefficients are assigned to.
%    A string specifying the type of geometric domain the coefficients
%    are assigned to. This string has two possible values: 'cell'
%    or 'face'. 
    RegionType;
end

properties
% RegionID - ID of the geometric regions the coefficients are assigned to.
%    For 2-D each RegionID satisfies 0 < RegionID(j) < NumFaces in the 
%    geometry. For 3-D each RegionID satisfies 0 < RegionID(j) < NumCells 
%    in the geometry.  
    RegionID;

% m - Second-order derivative coefficient
%     A numeric scalar, vector or matrix, or function
%     handle representing the m coefficient of the PDE.
     b;
     d;
     c;
     a;

end

methods (Hidden=true, Access={?pde.PDEModel})
    function obj=CoefficientAssignment(car, varargin) 
      obj.RecordOwner = car; 
      parser = inputParser;
      parser.PartialMatching=false; % Clash between 'face' and 'f'     
      parser.addParameter('face', []);     
      parser.addParameter('cell', []);     
      parser.addParameter('b', []);  
      parser.addParameter('d', []); 
      parser.addParameter('c', []); 
      parser.addParameter('a', []); 
      parser.addParameter('SystemSize', 1); 
      parser.parse(varargin{:});
           
      numdims = 2;
      if ~isempty(parser.Results.face)
              obj.RegionType = 'face';
              obj.RegionID = parser.Results.face;
      elseif ~isempty(parser.Results.cell)
             obj.RegionType = 'cell';
             obj.RegionID = parser.Results.cell; 
             numdims = 3;
      end       
      systemsize = parser.Results.SystemSize;

      obj.b = parser.Results.b;     
      obj.d = parser.Results.d; 
      obj.c = parser.Results.c; 
      obj.a = parser.Results.a; 
      obj.checkAllMatrixCoefSizes(systemsize, numdims);  
      obj.checkFcnHdlArgCounts(systemsize, numdims);
      obj.checkSparseness();
      % obj.checkNumericComplexity();
    end
end

methods
    
     function set.RegionType(self, rtype)     
%        self.ValidateRegionType(rtype);     
        self.RegionType = rtype;       
    end
    
    function set.RegionID(self, rids)      
        self.ValidateRegionID(rids);      
        self.RegionID = rids;        
    end        
    
    function set.b(self, coef)    
      self.CoefPrecheck(coef);     
      self.b = coef;     
    end
    
    function set.d(self, coef)    
      self.CoefPrecheck(coef);     
      self.d = coef;     
    end
    
    function set.c(self, coef)    
      self.CoefPrecheck(coef);     
      self.c = coef;     
    end
    
    function set.a(self, coef)    
      self.CoefPrecheck(coef);     
      self.a = coef;     
    end

end

methods(Hidden=true, Access = {?hhxt.CoefficientAssignmentRecords})
    function tf=sameCoefficients(self,other,varargin)
      if ~coefficientsMatch(self, other)
          tf = false;
          return   
      end
      if isempty(varargin)
             tf = (isequal(self.b, other.b) && ...
                   isequal(self.d, other.d) && ...
                   isequal(self.c, other.c) && ...
                   isequal(self.a, other.a) );                       
      else
          cname = varargin{1};
          switch cname
            case 'b'
                tf = isequal(self.b, other.b);  
            case 'd'
                tf = isequal(self.d, other.d);
            case 'c'
                tf = isequal(self.c, other.c); 
            case 'a'
                tf = isequal(self.a, other.a);
          end
      end
    end   
    
    function tf=numericCoefficients(self,varargin) 
      if isempty(varargin)
          
             tf = ~any(isa(self.b, 'function_handle')  || ...
                      isa(self.d, 'function_handle')  || ...
                      isa(self.c, 'function_handle')  || ...
                      isa(self.a, 'function_handle') );                                
      else
          cname = varargin{1};
          switch cname
            case 'b'
                tf = ~isa(self.b, 'function_handle');
            case 'd'
                tf = ~isa(self.d, 'function_handle');
            case 'c'
                tf = ~isa(self.c, 'function_handle');
            case 'a'
                tf = ~isa(self.a, 'function_handle');
          end
      end
    end 
            
    
    function tf = coefficientsMatch(self, other)
        tf = true;
        tf = tf & (self.bDefined() == other.bDefined());  
        tf = tf & (self.dDefined() == other.dDefined());         
        tf = tf & (self.cDefined() == other.cDefined());
%         tf = tf & (self.aDefined() == other.aDefined())
    end
    
    function performSolverPrecheck(self, systemsize, numfaces, numcells)
         ndims = 2;       
         if strcmp(self.RegionType, 'face')
            if any(self.RegionID > numfaces)
                error(message('pde:pdeCoefficientSpecification:invalidFaceIndexPresolve'));
            end
         else
            if any(self.RegionID > numcells)
                error(message('pde:pdeCoefficientSpecification:invalidCellIndexPresolve'));
            end 
            ndims = 3;
         end  
         checkAllMatrixCoefSizes(self, systemsize, ndims); 
         checkMandD(self);
         checkSparseness(self);
         checkFcnHdlArgCounts(self, systemsize, ndims)   
         % checkNumericComplexity(self); 
     end                
    
end


methods(Hidden=true, Access = private)        
    % checkAllMatrixCoefSizes - check the size of all coefficients that
    % are defined by a matrix.
    function checkAllMatrixCoefSizes(self, systemsize, ndims)
       self.checkBCoefSize(self.b, systemsize,ndims); 
       self.checkDCoefSize(self.d, systemsize);
       self.checkCCoefSize(self.c, systemsize,ndims);
       self.checkACoefSize(self.a, systemsize);
    end
    
    function checkSparseness(self)
        sparsecoef = issparse(self.b) | issparse(self.d) | issparse(self.c) | issparse(self.a);
        if sparsecoef
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueSparse'));   
        end        
    end  
    
    function checkFcnHdlArgCounts(self, systemsize, ndims)        
        if self.bDefined() 
            self.checkCoefFcnHdlArgCounts(self.b, systemsize, ndims);              
        end                           
        if self.dDefined() 
            self.checkCoefFcnHdlArgCounts(self.d, systemsize, ndims);         
        end               
        if self.cDefined()
            self.checkCoefFcnHdlArgCounts(self.c, systemsize, ndims);     
        end               
        if self.aDefined() 
            self.checkCoefFcnHdlArgCounts(self.a, systemsize, ndims);      
        end 
    end
    
end

methods(Hidden=true, Access = {?hhxt.CoefficientAssignmentRecords})
     function tf = bDefined(self)
        tf = self.coefDefined(self.b); 
     end
     
     function tf = dDefined(self)
        tf = self.coefDefined(self.d); 
     end

     function tf = cDefined(self)
        tf = self.coefDefined(self.c); 
     end

     function tf = aDefined(self)
        tf = self.coefDefined(self.a); 
     end
    
   function tf = hasComplexCoefficient(self, loc, state)  
        import hhxt.CoefficientAssignment.*
        tf = false(4,1);        
        tf(1) = coefIsComplexNumericOrFcnHdl(self.b, loc, state);   
        tf(2) = coefIsComplexNumericOrFcnHdl(self.d, loc, state);
        tf(3) = coefIsComplexNumericOrFcnHdl(self.c, loc, state);
        tf(4) = coefIsComplexNumericOrFcnHdl(self.a, loc, state);
   end    
   
end

methods(Static, Hidden=true)       
    % preDelete - Preserved for backward compatibility. In R2017b and prior
    % releases, the coefficient assignment has an ObjectBeingDestroyed 
    % callback that calls this function.
    function preDelete(self,~)
        if ~isempty(self.RecordOwner) && isvalid(self.RecordOwner)
            self.RecordOwner.delistCoefficientAssignment(self);
        end
    end    
end

methods(Hidden=true)
  function delete(self)
       if ~isempty(self.RecordOwner) && isvalid(self.RecordOwner)
            self.RecordOwner.delistCoefficientAssignment(self);
       end
  end
end

methods(Static, Access = private)
    function tf = coefDefined(coef)
       tf = (isnumeric(coef) && ~(isscalar(coef) && coef == 0) || isa(coef,'function_handle'));
    end

    function tf = coefIsComplexNumericOrFcnHdl(coef, loc, state)
         tf = false;
         if ~hhxt.CoefficientAssignment.coefDefined(coef)
             return 
         end
         if isnumeric(coef)
            if ~isreal(coef)
                tf = true;          
            end 
         else % isa(coef, 'function_handle')
             res = coef(loc, state);
             if ~isreal(res)
                tf = true;          
             end 
         end              
    end        

    function ok=ValidateRegionID(rgnid)
      % Must be real(non-complex), full, natural number.        
      if ~isreal(rgnid) || ~all(rgnid(:) > 0) || issparse(rgnid) || any(mod(rgnid(:),1)~=0)
        error(message('pde:pdeCoefficientSpecification:invalidRegionID'));   
      end      
      ok = true; 
    end
    function ok=CoefPrecheck(coef)
      if isfloat(coef)                   
        if any(isnan(coef))
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueNaN'));
        elseif any(isinf(coef))
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueInf'));
        elseif(isempty(coef))
            error(message('pde:pdeCoefficientSpecification:invalidCoefValueEmpty'));
        end
      elseif ~isa(coef, 'function_handle')            
          if pde.CoefficientAssignment.stringCoefMissedDef(coef)
            error(message('pde:pdeCoefficientSpecification:missedCoefValue'));   
          elseif ischar(coef)
             error(message('pde:pdeCoefficientSpecification:invalidCoefValueString')); 
          else
            error(message('pde:pdeCoefficientSpecification:invalidCoefValue'));       
          end                               
      end      
      ok = true; 
    end
    
    function checkBCoefSize(bcoef, systemsize, ndims)
        if isscalar(bcoef) && ~isa(bcoef, 'function_handle') && bcoef == 0
            return
        elseif isa(bcoef, 'function_handle')
            return
        end
        cveclen = numel(bcoef);      
        if ndims == 2
             lengthok = (cveclen == 2*systemsize);
             if ~(isvector(bcoef)  && iscolumn(bcoef) && lengthok)      
                error(message('pde:pdeCoefficientSpecification:invalidCMatrixSize'));   
             end                    
        else
            lengthok = (cveclen == 3*systemsize);
             if ~(isvector(bcoef)  && iscolumn(bcoef) && lengthok)      
                error(message('pde:pdeCoefficientSpecification:invalidCMatrixSize'));   
             end               
        end                      
    end
    
    function checkDCoefSize(dcoef, systemsize)
         if isscalar(dcoef) && ~isa(dcoef, 'function_handle') && dcoef == 0
            return
         elseif isa(dcoef, 'function_handle')
             return
         end
         dveclen = numel(dcoef);
         lengthok = (dveclen == 0 || dveclen == 1 || dveclen == systemsize || ...
                     dveclen == (systemsize*(systemsize+1)/2) || dveclen == systemsize*systemsize);          
         if ~(isvector(dcoef)  && iscolumn(dcoef) && lengthok)      
            error(message('pde:pdeCoefficientSpecification:invalidDMatrixSize'));   
         end
    end
    function checkCCoefSize(ccoef, systemsize, ndims)
        if isscalar(ccoef) && ~isa(ccoef, 'function_handle') && ccoef == 0
            return
        elseif isa(ccoef, 'function_handle')
            return
        end
        cveclen = numel(ccoef);      
        if ndims == 2
             lengthok = (cveclen == 1 || cveclen == 2 || cveclen == 3 || ...
                         cveclen == 4 || cveclen == systemsize || ...
                         cveclen == 2*systemsize || cveclen == 3*systemsize || ...
                         cveclen == 4*systemsize || cveclen == 2*systemsize*((2*systemsize)+1)/2 || ...
                         cveclen == 4*(systemsize^2));
             if ~(isvector(ccoef)  && iscolumn(ccoef) && lengthok)      
                error(message('pde:pdeCoefficientSpecification:invalidCMatrixSize'));   
             end                    
        else
            lengthok = (cveclen == 1 || cveclen == 3 || cveclen == 6 || ...
                         cveclen == 9 || cveclen == systemsize || ...
                         cveclen == 3*systemsize || cveclen == 6*systemsize || ...
                         cveclen == 9*systemsize || cveclen == 3*systemsize*((3*systemsize)+1)/2 || ...
                         cveclen == 9*(systemsize^2));
             if ~(isvector(ccoef)  && iscolumn(ccoef) && lengthok)      
                error(message('pde:pdeCoefficientSpecification:invalidCMatrixSize'));   
             end               
        end                      
    end
    function checkACoefSize(acoef, systemsize) 
         if isscalar(acoef) && ~isa(acoef, 'function_handle') && acoef == 0
            return
         elseif isa(acoef, 'function_handle') 
            return
         end
         aveclen = numel(acoef);        
         lengthok = (aveclen == 0 || aveclen == 1 || aveclen == systemsize || ...
                     aveclen == (systemsize*(systemsize+1)/2) || aveclen == systemsize*systemsize);          
         if ~(isvector(acoef)  && iscolumn(acoef) && lengthok)      
            error(message('pde:pdeCoefficientSpecification:invalidAMatrixSize'));   
        end
    end  
    
    %
    % stringCoefMissedDef
    % Returns true if the coefficient specification is a char of length 1
    % and is either 'm', 'd', 'c','a', 'f'.
    % This detects a typo that's not uncommon;  e.g. 
    % (..., 'm', 1, 'd', 2, 'c', 'a', 'f', 0)
    function tf = stringCoefMissedDef(coef)
        tf = false;
        if ischar(coef) && numel(coef) == 1
            tf = strcmp('b',coef) || strcmp('d',coef) || strcmp('c',coef) ...
                 || strcmp('a',coef);
        end
    end           
         
    
%
%   Coefficient function-handle argument checks
%   
    function checkCoefFcnHdlArgCounts(coef, systemsize, ndims)
        if ~isa(coef, 'function_handle')
            return
        end
        location.x = 0;
        location.y = 0;          
        location.z = 0;
        location.subdomain=1;  
        state.u = zeros(systemsize, 1);
        state.ux = zeros(systemsize, 1);
        state.uy = zeros(systemsize, 1);
        state.uz = zeros(systemsize, 1);
        state.time = 0;  
        str = round((systemsize-1)/2);
        state.kx = zeros(str, 1);
        state.ky = zeros(str, 1);
        state.kz = zeros(str, 1);
        state.uDx = zeros(str, 1);
        state.uDy = zeros(str, 1);
        state.uDz = zeros(str, 1);
        state.Rex = zeros(str, 1);
        state.Rey = zeros(str, 1);
        state.Rez = zeros(str, 1);
        state.uDDELTx = zeros(str, 1);
        state.uDDELTy = zeros(str, 1);
        state.uDDELTz = zeros(str, 1);
        try 
            coef(location, state);
        catch
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));   
        end
        
        try 
            res = coef(location, state);
        catch
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));   
        end
              
        if ~isa(res,'double')
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleOutputType'));
        end
        
        if ndims == 2
           p = [0 1 0; 0 0 1];
           t = [1; 2; 3];          
        else            
           p = [0 1 0 0; 0 0 1 0; 0 0 0 1]; 
           t = [1; 2; 3; 4];           
        end
        u = 0;
        time = 0;
        threwerr = false; 
        try 
            res = coef(p,t,u,time);
        catch
            threwerr = true;
        end        
        if ~threwerr
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));    
        end

        if ~isa(res,'double')
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleOutputType'));
        end        
    end  
    
end    
   properties (Hidden = true, Access='private')
      RecordOwner;
   end   
end
