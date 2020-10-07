function ca = findCoefficients(self, varargin)  
% findCoefficients Find coefficients assignment for a geometric region
%    CA = findCoefficients(CC,REGIONTYPE, REGIONID) returns the 
%    CoefficientAssignment object CA that defines the equation coefficients
%    assigned to the given REGIONTYPE and REGIONID. Where REGIONTYPE is  
%    'cell' for a 3-D model or 'face' for a 2-D model and REGIONID is the 
%    ID of the geometric entity.
%  
%    See also pde.PDEModel/specifyCoefficients, NumericCoefficientFormat, 
%              FunctionCoefficientFormat
% 

% Copyright 2015-2016 The MathWorks, Inc.

narginchk(3,3);  
% if numassigns == 0
%     return;
% end

if isempty(self.ParentPdemodel.Geometry)
    return;
end

parser = inputParser;
parser.addParameter('cell', [], @pde.PDEModel.isValidEntityID); 
parser.addParameter('face', [], @pde.PDEModel.isValidEntityID);  
parser.parse(varargin{:});

if ~ismember('face', parser.UsingDefaults)
    if ~self.ParentPdemodel.IsTwoD
      error(message('pde:pdeCoefficientSpecification:noBoundaryCoefficients')); 
    end
    qregiontype = 'face';
    qregionid = parser.Results.face;
    maxid = self.ParentPdemodel.Geometry.NumFaces;       
    if any(qregionid > self.ParentPdemodel.Geometry.NumFaces)
       error(message('pde:pdeCoefficientSpecification:invalidRegionID'));
    end 
else
    qregiontype = 'cell';
    qregionid = parser.Results.cell;
    if self.ParentPdemodel.IsTwoD
        maxid = 0;
    else
        maxid = self.ParentPdemodel.Geometry.NumCells;
    end       
end

if any(qregionid > maxid)
       error(message('pde:pdeCoefficientSpecification:invalidRegionID'));
end    

numassigns = numel(self.CoefficientAssignments);
numqueries = numel(qregionid);
ca = repmat(self.CoefficientAssignments(1),size(qregionid));
for i = 1:numqueries
    rid = qregionid(i);
    cafound = false;
    for j = numassigns:-1:1
       thisca =  self.CoefficientAssignments(j);       
       if strcmp(qregiontype,thisca.RegionType) && ismember(rid, thisca.RegionID)        
          ca(i) = thisca;
          cafound = true;
          break
       end
    end
    if ~cafound
      error(message('pde:pdeCoefficientSpecification:entitiesWithoutCoefficients'));
    end   
end

end

% LocalWords:  REGIONTYPE REGIONID
