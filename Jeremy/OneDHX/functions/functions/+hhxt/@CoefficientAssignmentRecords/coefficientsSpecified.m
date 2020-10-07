function tf = coefficientsSpecified(self, varargin)  
% coefficientsSpecified true if coefficients are specified for domain or subdomain
%    TF = coefficientsSpecified(CR,'face', FACEID) returns true if coefficients 
%    are specified for a face whos ID is given by FACEID.
%  
%    TF = coefficientsSpecified(CR) returns true if coefficients are specified 
%    on all subdomains in the domain.
%  
%    See also pde.PDEModel/specifyCoefficients, NumericCoefficientFormat, 
%              FunctionCoefficientFormat
% 

% Copyright 2015-2016 The MathWorks, Inc.


narginchk(1,3);  
tf = false;
parser = inputParser;
parser.addParameter('face', [], @pde.PDEModel.isValidEntityID);   
parser.parse(varargin{:});

if isempty(self.ParentPdemodel.Geometry)
    return;
end

if isempty(parser.UsingDefaults) 
    regionid = parser.Results.face;
    if any(regionid > self.ParentPdemodel.Geometry.NumFaces)
       error(message('pde:pdeCoefficientSpecification:invalidRegionID'));
    end      
    tf = false(size(regionid));
end




if ~(self.ParentPdemodel.IsTwoD)
    tf = true; % Can only be global for 3-D
    return;
end

globalquery=false;
% Face query
if ~isempty(parser.UsingDefaults) && strcmp(parser.UsingDefaults{1},'face')
    regionid = 1:self.ParentPdemodel.Geometry.NumFaces;
    tf = false(size(regionid));
    globalquery=true;
end

numqueries = numel(regionid);
numassigns = numel(self.CoefficientAssignments);
for i = 1:numqueries
    rid = regionid(i);
    for j = 1:numassigns
       thisca =  self.CoefficientAssignments(j);
       if ismember(rid, thisca.RegionID) 
          tf(i) = true;
          break       
       end
    end    
end

if globalquery
   tf = all(tf == true);    
end

end

% LocalWords:  subdomain FACEID subdomains
