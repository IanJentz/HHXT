function tf = globalCoefficients(self, varargin)  
% globalCoefficients  True if consistent across the geometric domain
%    tf = globalCoefficients(CC) returns TRUE if all coefficients are 
%    consistent across the geometric domain, otherwise returns FALSE. 
%  
%    tf = globalCoefficients(CC, COEFNAME) returns TRUE if the coefficient of
%    type COEFNAME is consistent across the geometric domain, otherwise 
%    returns FALSE.
%  
%   See also  pde.PDEModel/specifyCoefficients

% Copyright 2015-2017 The MathWorks, Inc.

narginchk(1,2); 
tf = true;
% numassigns = numel(self.CoefficientAssignments);
% if numassigns == 0
%     tf = false;
%     return;
% end

if ~self.coefficientsSpecified()   
    tf = false;
    return;
end

if self.ParentPdemodel.Geometry.NumCells==1
    % 3-D geometry with one cell
    return;
elseif self.ParentPdemodel.Geometry.NumFaces==1
    % 2-D geometry with one face
    return;
end

if self.ParentPdemodel.IsTwoD
    qregiontype = 'face';
    nRegion = self.ParentPdemodel.Geometry.NumFaces;
else
    qregiontype = 'cell';
    nRegion = self.ParentPdemodel.Geometry.NumCells;
end



% if isempty(self.ParentPdemodel.Geometry)
%     % Coefficients defined, but no geometry
%     return;
% end



if nargin == 2
   coefname =  varargin{1};
   iscoef = false;
   if (ischar(coefname) && numel(coefname)==1)
       iscoef = (strcmp(coefname, 'b'));
   end
   if ~iscoef
      error(message('pde:pdeCoefficientSpecification:invalidCoefName')); 
   end
end
    

systemsize = self.ParentPdemodel.PDESystemSize;
location.x = 0;
location.y = 0;
location.z = 0;
location.subdomain = 1;
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
    
thiscoef = self.findCoefficients(qregiontype, 1);
thisBVal = evaluateFcnHdlVal(thiscoef.b,location,state);

for nsd = 2:nRegion
    othercoef = self.findCoefficients(qregiontype, nsd);           
    
if (~self.allCoefficientsNumeric)
            location.subdomain = nsd;
            otherBVal = evaluateFcnHdlVal(othercoef.b,location,state);
            if (...
                    ( any(thisBVal(:) ~= otherBVal(:))) ...
                )
                tf = false;
                return;
            end
    
else
   
    
    if ~thiscoef.sameCoefficients(othercoef,varargin{:})
        tf = false;
        return;
    end            
end  
end



end


% Helper function to find values of coefficients returned from user defined
% function.
  function res = evaluateFcnHdlVal(coef, location, state)
        if ~isa(coef, 'function_handle')
            res = 0;
            return
        end
        
        try
            res = coef(location, state);
        catch
            error(message('pde:pdeCoefficientSpecification:invalidFcnHandleArgs'));
        end
  end    
    
  