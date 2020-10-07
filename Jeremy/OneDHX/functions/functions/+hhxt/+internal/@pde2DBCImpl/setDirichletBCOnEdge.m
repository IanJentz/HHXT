function [h, r] = setDirichletBCOnEdge(self, bc, femEdges, h, r)
%setDirichletBCOnEdge Evaluate h and r for Dirichlet BC
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.
numFemEdges = length(femEdges);
ne = self.numEdges;
% Handle the first term (r vector)
if(~isempty(bc.term1))
  if isa(bc.term1, 'function_handle')
    [riV1, riV2] = self.callDirichletFuncOnEdge(bc, femEdges, bc.term1, 1);
    r(:,femEdges) = riV1; r(:,femEdges+ne) = riV2;
  else
    ri = bc.term1(:);
    % For Dirichlet, the length of term1 must equal N
    if(length(ri) ~= self.N)
      error(message('pde:pde2DBCImpl:invalidLenR', bc.appRegionID, self.N));
    end
    rAllEdges = repmat(ri, 1, numFemEdges);
    r(:,femEdges) = rAllEdges; r(:,femEdges+ne) = rAllEdges;
  end
end
% Handle the second term (h matrix)
if(isempty(bc.term2))
  hAllEdges = repmat(self.hDefault, 1, numFemEdges);
  h(:,femEdges) = hAllEdges; h(:,femEdges+ne) = hAllEdges;
else
  if isa(bc.term2, 'function_handle')
    [hV1, hV2] = self.callDirichletFuncOnEdge(bc, femEdges, bc.term2, 2);
    h(:,femEdges) = hV1; h(:,femEdges+ne) = hV2;
  else
    hi = bc.term2(:);
    % For Dirichlet, the length of term2 must equal N^2
    if(length(hi) ~= self.N^2)
      error(message('pde:pde2DBCImpl:invalidLenH', bc.appRegionID, self.N^2));
    end
    hAllEdges = repmat(hi, 1, numFemEdges);
    h(:,femEdges) = hAllEdges; h(:,femEdges+ne) = hAllEdges;
  end
end

end

