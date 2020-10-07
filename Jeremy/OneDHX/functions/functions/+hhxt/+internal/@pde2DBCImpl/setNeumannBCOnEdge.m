function [q, g] = setNeumannBCOnEdge(self, bc, femEdges, q, g)
%setNeumannBCOnEdge Evaluate q and g coefficients for a Neumann BC
% The input may be either numeric or a user-written function.
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.
numFemEdges = length(femEdges);

% Handle the first term
if(~isempty(bc.term1))
  if isa(bc.term1, 'function_handle')
    gi = self.callNeumannFuncOnEdge(bc, femEdges, bc.term1, 1);
    g(:,femEdges) = gi;
  else
    if(length(bc.term1) ~= self.N)
      error(message('pde:pde2DBCImpl:invalidLenG', bc.appRegionID, self.N));
    end
    t1 = bc.term1(:);
    gAllEdges = repmat(t1, 1, numFemEdges);
    g(:,femEdges) = gAllEdges;
  end
end

% Handle the second term
if(~isempty(bc.term2))
  if isa(bc.term2, 'function_handle')
    qi = self.callNeumannFuncOnEdge(bc, femEdges, bc.term2, 2);
    q(:,femEdges) = reshape(qi, [], numFemEdges);
  else
    if(length(bc.term2(:)) ~= self.N^2)
      error(message('pde:pde2DBCImpl:invalidLenQ', bc.appRegionID, self.N^2));
    end
    qAllEdges = repmat(bc.term2(:), 1, numFemEdges);
    q(:,femEdges) = qAllEdges;
  end
end
end

