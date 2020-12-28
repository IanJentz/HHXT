function self = AdaptiveMeshRefinement(self,R)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


elem_refine = [];

nodes = find(any([...
    R.ConvAbsolute >= self.SolverOptions.AbsoluteTolerance,...
    R.ConvRelative >= self.SolverOptions.RelativeTolerance,...
    R.ConvResidual >= self.SolverOptions.ResidualTolerance,...
    ],2));
if ~isempty(nodes)
elem_N = self.Mesh.findElements('attached',nodes);
elem_refine = [elem_refine;elem_N']; % need to be passed as a column
end

% if isempty(R.ZGradients)
%     MGradient = sqrt(R.XGradients.^2+R.YGradients.^2);
% else
%     MGradient = sqrt(R.XGradients.^2+R.YGradients.^2+R.ZGradients.^2);
% end
% 
% max_refine_nodes = floor(0.05*size(MGradient,1));
% 
% for j = 1:n1
%     [~,nodes] = sort(MGradient(:,j),'descend');
%     nodes = nodes(1:max_refine_nodes);
%     elem_N = self.Mesh.findElements('attached',nodes);
%     elem_refine = [elem_refine;elem_N']; % need to be passed as a column
% end

for j = 1:size(R.ElemHeating,2)
    nodes = find(...
        abs(R.ElemHeating(:,j)-mean(R.ElemHeating(:,j)))...
        >= 2*std(R.ElemHeating(:,j)) );
    elem_N = self.Mesh.findElements('attached',nodes);
    elem_refine = [elem_refine;elem_N']; % need to be passed as a column
end

if ~isempty(elem_refine)
um = R.NodalSolution;%reshape(u,size(u,1)/n1,n1);
[self,u0m] = refineHHXTmesh(self,um,elem_refine);  % elements to refine, input as a column
u0 = u0m(:);

R.NodalSolution = u0;
setInitialConditions(self,R);

% %prebuild the Gauss Quadrature points
% self = self.modelQuadPoints();
% % prepare models non linear variables
% self = updateQuadPointVals(self,u0); % pre-build gauss quad points
% self = self.updatePerm;  % update the permeability at all gauss points
end


end

