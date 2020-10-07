function tf = boundaryConditionsSpecified(self, varargin)  
%boundaryConditionsSpecified returns true if BCs are specified for boundaries
%   TF = boundaryConditionsSpecified(BCR,REGIONTYPE, REGIONID) returns true if 
%   Boundary conditions are specified on a region of the doamin defined by 
%   REGIONTYPE and REGIONID. Where REGIONTYPE is one of the following: 
%   'face' or 'edge' and REGIONID is the corresponding ID.
%  
%   See also pde.PDEModel/applyBoundaryCondition
% 
% Copyright 2016 The MathWorks, Inc.

narginchk(3,3);  

if isempty(self.ParentPdemodel.Geometry)
    return;
end

parser = inputParser;
parser.addParameter('face', [], @pde.PDEModel.isValidEntityID); 
parser.addParameter('edge', [], @pde.PDEModel.isValidEntityID);  
parser.parse(varargin{:});

if ~ismember('face', parser.UsingDefaults)
    qregiontype = 'face';
    qregionid = parser.Results.face;
    maxid = self.ParentPdemodel.Geometry.NumFaces;       
    if any(qregionid > maxid)
       error(message('pde:pdeModel:invalidFaceIndex'));
    end 
    if self.ParentPdemodel.IsTwoD
      error(message('pde:pdeBoundaryConditions:noBoundaryCoefficients')); 
    end
else
    qregiontype = 'edge';
    qregionid = parser.Results.edge;
    maxid = self.ParentPdemodel.Geometry.NumEdges;
    if any(qregionid > maxid)
       error(message('pde:pdeModel:invalidEdgeIndex'));
    end 
    if ~self.ParentPdemodel.IsTwoD
      error(message('pde:pdeBoundaryConditions:noBoundaryCoefficients')); 
    end
end

numassigns = numel(self.BoundaryConditionAssignments);
numqueries = numel(qregionid);
tf = false(size(qregionid)); 
for i = 1:numqueries
    rid = qregionid(i);
    for j = numassigns:-1:1
       thisbc =  self.BoundaryConditionAssignments(j);       
       if strcmpi(qregiontype,thisbc.RegionType) && ismember(rid, thisbc.RegionID)       
          tf(i) = true;
          break
       end
    end
end

end