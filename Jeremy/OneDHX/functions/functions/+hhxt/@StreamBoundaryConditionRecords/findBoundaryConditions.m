function bc = findBoundaryConditions(self, varargin)  
% findBoundaryConditions finds boundary condition assignment for a boundary
%    BC = findBoundaryConditions(BCR,REGIONTYPE,REGIONID) returns the 
%    BoundaryCondition object BC that defines the boundary conditions
%    assigned to the given REGIONTYPE and REGIONID. Here REGIONTYPE is  
%    'Face' for a 3-D model or 'Edge' for a 2-D model and REGIONID is the 
%    ID of the geometric entity.
%  
%    See also pde.PDEModel/applyBoundaryCondition
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
bc = repmat(self.BoundaryConditionAssignments(1),size(qregionid));
for i = 1:numqueries
    rid = qregionid(i);
    bcfound = false;
    for j = numassigns:-1:1
       thisbc =  self.BoundaryConditionAssignments(j);       
       if strcmpi(qregiontype,thisbc.RegionType) && ismember(rid, thisbc.RegionID)       
          bc(i) = thisbc;
          bcfound = true;
          break
       end
    end
    if ~bcfound
        error(message('pde:pdeBoundaryConditions:entitiesWithoutBCs'));
    end   
end

end