classdef(Hidden) pde2DBCSfnImpl < hhxt.internal.pdeBCImpl
    %pde2DBCSfnImpl class for 2-D linear and quadratic element boundary conditions
    %
    % This undocumented function may be changed or removed in a future release.
    
    %       Copyright 2016 The MathWorks, Inc.
    
    properties(Access=private)
        edgeBoundaryConditions;
        points,  numPoints, time, uN, gotTime, gotU;
    end
    
    methods
        function obj = pde2DBCSfnImpl(myPde,bc,isR2016b)
            obj@pde.internal.pdeBCImpl(myPde);
            obj.copyBoundaryConditionsArray(bc,isR2016b);
        end
        
        function [Q,G,H,R] = getBCMatrices(self,p,e,u,time,meshorder)
            %getBCMatrices Calculate global matrices for linear and
            %quadratic 2-D boundary conditions
            %
            % This function iterates over all BCs applied to the edge of a
            % 2-D model. First it gets the element edge (also referred to
            % as discrete edge) for the particular geometry edge. Then it
            % calls the appropriate lower-level function to compute the BC
            % terms for each element edge.
            
            import pde.internal.*;
            
            self.points = p;
            self.numPoints = size(p,2);
            neqn = self.N*self.numPoints;
            
            % Get rid of unwanted edges, that are part of internal
            % boundaries between subdomians.
            ie=pdesde(e);
            e=e(:,ie);
            
            self.gotU = nargin > 3 && ~isempty(u);
            self.gotTime = nargin > 4 && ~isempty(time);
            
            if(self.gotU)
                self.uN = reshape(u,[],self.N)';
            else
                self.uN = zeros(self.N, self.numPoints);
            end
            
            if(self.gotTime)
                self.time = time;
            else
                self.time = [];
            end
            
            Q = sparse(neqn,neqn);
            G = sparse(neqn,1);
            
            % map of Dirichlet BCs with node ID as the key
            dirBoundaryConditions = containers.Map('KeyType','uint32','ValueType','any');
            
            % map of neumann BCs with geom edge ID as key.
            % This "map" just tracks which edges already have Neumann BCs applied to
            % them so they are applied only once.
            neuBoundaryConditions = containers.Map('KeyType','uint32','ValueType','logical');
            
            % iterate, last to first so only last Neumann BC on a edge is applied
            for i=length(self.edgeBoundaryConditions):-1:1
                bcsi = self.edgeBoundaryConditions(i);
                %each column bcsi.elemEdges will contains node IDs of the
                %FE edges associated with feature edge given by
                %bcsi.appRegionID.
                if strcmp(meshorder,'quadratic')
                    bcsi.elemEdges = extractDiscreteEdgeQuadratic(self, e, bcsi.appRegionID);
                else
                    bcsi.elemEdges = extractDiscreteEdgeLinear(self, e, bcsi.appRegionID);
                end
                
                % Call the appropriate BC function for each type
                if(bcsi.bcType == pdeBCImpl.valueBC)
                    bcsi.nodes = self.getNodesForEdges(bcsi.elemEdges);
                    bcsi.r = bcsi.term1;
                    bcsi.h = bcsi.term2;
                    self.setValueBCOnEdge(bcsi, dirBoundaryConditions);
                elseif(bcsi.bcType == pdeBCImpl.neumannBC)
                    % Last Neumann BC on a edge is the right one. So, insert only the last one.
                    if ~neuBoundaryConditions.isKey(bcsi.edgeID)
                        neuBoundaryConditions(bcsi.edgeID) = true;
                        bcsi.g = bcsi.term1;
                        bcsi.q = bcsi.term2;
                        [Qi, Gi] = setNeumannBCOnEdge(self, bcsi);
                        Q = Q + Qi;
                        G = G + Gi;
                    end
                elseif(bcsi.bcType == pdeBCImpl.dirichletBC)
                    bcsi.nodes = self.getNodesForEdges(bcsi.elemEdges);
                    bcsi.r = bcsi.term1;
                    bcsi.h = bcsi.term2;
                    setDirichletBCOnEdge(self, bcsi, dirBoundaryConditions);
                else
                    error(message('pde:pde2DBCImpl:invalidBCType', bcsi.bcType));
                end
            end
            
            %
            % Create the H and R matrices
            %
            % preallocate some space to hold the triplets
            numDirBCs = length(dirBoundaryConditions);
            preAllocSize = max(numDirBCs*self.N, 1);
            tripR = zeros(2, preAllocSize);
            tripH = zeros(3, preAllocSize);
            bcNodes = dirBoundaryConditions.keys;
            nodalVals = dirBoundaryConditions.values;
            % actual counts of constraints and terms
            consCount = 0; hCount = 0; rCount = 0;
            for i=1:length(dirBoundaryConditions)
                nodeID = bcNodes{i};
                hr = nodalVals{i};
                hNodeI = hr.h;
                rNodeI = hr.r;
                for r=1:self.N
                    if(any(hNodeI(r,:)~=0))
                        rr = rNodeI(r);
                        consCount = consCount + 1;
                        if(rr~=0)
                            rCount = rCount + 1;
                            tripR(1, rCount) = consCount;
                            tripR(2, rCount) = rr;
                        end
                        for c=1:self.N
                            hrc = hNodeI(r,c);
                            if(hrc~=0)
                                eq = (c-1)*self.numPoints + nodeID;
                                hCount = hCount + 1;
                                tripH(1, hCount) = consCount;
                                tripH(2, hCount) = eq;
                                tripH(3, hCount) = hrc;
                                hNodeI(r,c) = hrc;
                            end
                        end
                    end
                end
            end
            
            
            if(consCount)
                H = sparse(tripH(1, 1:hCount), tripH(2, 1:hCount), tripH(3, 1:hCount), consCount, neqn);
                R = sparse(tripR(1, 1:rCount), 1, tripR(2, 1:rCount), consCount, 1);
            else
                H = sparse(1,neqn);
                R = sparse(neqn,1);
            end
            
        end
        
    end % methods
    
    methods(Access=private)
        function self = copyBoundaryConditionsArray(self,bc,isR2016b)
            % copy pdeBoundaryConditions to local data structure
            
            import pde.internal.pdeBCImpl;
            
            if(~isempty(self.myPde))
                self.edgeBoundaryConditions = struct('bcType',{},'term1',{},'term2',{}, ...
                    'isVectorized',{},'appRegionID',{}, 'edgeID',{});
                if isempty(bc)
                    self.edgeBoundaryConditions = [];
                elseif isR2016b
                    for i=1:length(bc)
                        bci = bc(i);
                        val = struct();
                        if(strcmpi(bci.Type,'neumann'))
                            val.bcType = pdeBCImpl.neumannBC;
                            val.term1 = bci.g;
                            val.term2 = bci.q;
                        elseif(strcmpi(bci.Type,'dirichlet'))
                            val.bcType = pdeBCImpl.dirichletBC;
                            val.term1 = bci.r;
                            val.term2 = bci.h;
                        elseif(strcmpi(bci.Type,'value'))
                            val.bcType = pdeBCImpl.valueBC;
                            val.term1 = bci.u;
                            val.term2 = bci.EquationIndex;
                        elseif(strcmpi(bci.Type,'mixedvalue'))
                            val(1).bcType = pdeBCImpl.valueBC;
                            val(1).term1 = bci.u;
                            val(1).term2 = bci.EquationIndex;
                            val(2).bcType = pdeBCImpl.neumannBC;
                            val(2).term1 = bci.g;
                            val(2).term2 = bci.q;
                        elseif(strcmpi(bci.Type,'mixeddirichlet'))
                            val(1).bcType = pdeBCImpl.dirichletBC;
                            val(1).term1 = bci.r;
                            val(1).term2 = bci.h;
                            val(2).bcType = pdeBCImpl.neumannBC;
                            val(2).term1 = bci.g;
                            val(2).term2 = bci.q;
                        else
                            error(message('pde:pde2DBCImpl:invalidBCType', bci.Type));
                        end
                        for j = 1:numel(val)
                            val(j).isVectorized = strcmpi(bci.Vectorized, 'on');
                            val(j).appRegionID = i;
                            val(j).edgeID = i;
                            self.edgeBoundaryConditions(end+1) = val(j);
                        end
                    end
                else
                    for i=1:length(bc)
                        bci = bc(i);
                        % A case with no parameters is a zero-Neumann BC
                        % This is the default for FEM so we just skip it
                        if isempty(bci.Type)
                            continue;
                        end
                        if(strcmpi(bci.Type,'neumann'))
                            val.bcType = pdeBCImpl.neumannBC;
                            val.term1 = bci.g;
                            val.term2 = bci.q;
                        elseif(strcmpi(bci.Type, 'dirichlet')||strcmpi(bci.Type, 'mixeddirichlet'))
                            val.bcType = pdeBCImpl.dirichletBC;
                            val.term1 = bci.r;
                            val.term2 = bci.h;
                        elseif(strcmpi(bci.Type, 'value')||strcmpi(bci.Type, 'mixedvalue'))
                            val.bcType = pdeBCImpl.valueBC;
                            val.term1 = bci.u;
                            val.term2 = bci.EquationIndex;
                        else
                            error(message('pde:pde2DBCImpl:invalidBCType', bci.Type));
                        end
                        val.isVectorized = strcmpi(bci.Vectorized, 'on');
                
                        % A pdeBoundaryCondition can have multiple faces in its application region
                        % Separate them out here.
                        
                        for j=1:numel(bci.RegionID)
                            val.appRegionID = bci.RegionID(j);
                            if strcmp(bci.RegionType, 'Edge')
                                jj = bci.RegionID(j);
                                val.edgeID = jj;
                                self.edgeBoundaryConditions(end+1) = val;
                            end
                        end
                    end
                end
            end
        end
        
        
        function setValueBCOnEdge( self, bci, dirBoundaryConditions )
            %setValueBCOnFedge Calculate r and h entries for a value type BC
            % This function handles the case when a 'u' and/or 'EquationIndex' parameter
            % has been included in a pdeBoundaryConditions entity.
            
            import pde.internal.*;
            
            edgeID = bci.appRegionID;
            numNodes = length(bci.nodes);
            
            % Both the first and second terms are optional. But to get to this point,
            % at least one must have been included to determine the type.
            
            nodesCell = num2cell(bci.nodes);
            doesNodeHaveBC = dirBoundaryConditions.isKey(nodesCell);
            
            % Following flags are used to check if any node on this FE edge
            % or any DoF at a node still needs BC assignment. If not, do
            % nothing and return to the caller.
            
            nodesNeedBCAssignment = false;
            dofsNeedBCAssignment = false;
            
            if (~all(doesNodeHaveBC))
                nodesNeedBCAssignment = true;
            end
            
            if (nodesNeedBCAssignment == false)
                % check for BC at all DoFs only when all nodes already have some BC
                for n = cell2mat(nodesCell(doesNodeHaveBC))'
                    hr = dirBoundaryConditions(n);
                    previousBCAssignment = hr.r ~=0;
                    for j = 1:numel(hr.r)
                        if (~previousBCAssignment(j)) % If there was no previous assignment
                            dofsNeedBCAssignment = true;
                            continue
                        end
                    end
                end
            end
            
            if ~(nodesNeedBCAssignment || dofsNeedBCAssignment)
                return
            end
            
            isFuncU = isa(bci.term1, 'function_handle');
            
            % Check to see if a second term has been included.
            if(~isempty(bci.term2))
                % we have a second term
                if ~isnumeric(bci.term2)
                    error(message('pde:pde2DBCImpl:term2NotNum', edgeID));
                end
                t2 = bci.term2;
                rr1 = t2(:);
                lenRr1 = length(rr1);
                ri = zeros(self.N,1);
                rnz = ri~=0;
                if(lenRr1 > self.N)
                    eId = 'pde:pde2DBCImpl:t2TooLarge';
                    err = message(eId, lenRr1, self.N, edgeID);
                    throwAsCaller(MException(eId, err.getString()));
                end
                hMat = getH(self.N, t2, edgeID);
            else
                % no second term, default h-matrix (all dofs)
                rr1 = 1:self.N;
                lenRr1 = self.N;
                hMat = eye(self.N);
            end
            
            %
            % Now handle the first term
            %
            
            if isFuncU
                %if ~all(doesNodeHaveBC)
                edgeNormals = self.calcEdgeNormalsAtNodes(bci.elemEdges, bci.nodes);
                uVec = self.callValueFuncOnEdge(bci, bci.nodes, bci.term1, edgeNormals);
                nr = size(uVec, 1);
                if(nr > self.N)
                    eId = 'pde:pde2DBCImpl:t1TooLarge';
                    err = message(eId, nr, self.N, edgeID);
                    throwAsCaller(MException(eId, err.getString()));
                end
                checkTermLengths(nr, lenRr1, edgeID);
                rMat = zeros(self.N, numNodes);
                rMat(rr1,:) = uVec;
                %end
            elseif isnumeric(bci.term1)
                rMat = zeros(self.N, 1);
                if ~isempty(bci.term1)
                    lenT1 = length(bci.term1);
                    if(lenT1 > self.N)
                        eId = 'pde:pde2DBCImpl:t1TooLarge';
                        err = message(eId, lenT1, self.N, edgeID);
                        throwAsCaller(MException(eId, err.getString()));
                    end
                    checkTermLengths(lenT1, lenRr1, edgeID);
                    rMat(rr1) = bci.term1(:);
                    ri = rMat;
                    rnz = ri~=0;
                end
            else
                error(message('pde:pde2DBCImpl:term1NotNum', edgeID));
            end
            
            hrZ = HRNode(self.N);
            
            hi = hMat;
            hnz = hi~=0;
            
            
            for i=1:numNodes
                nodeID = bci.nodes(i);
                if(doesNodeHaveBC(i))
                    hr = dirBoundaryConditions(nodeID);
                else
                    hr = hrZ;
                end
                if(isFuncU)
                    ri = rMat(:,i);
                    rnz = ri~=0;
                end
                hr.nodeID = nodeID;
                hr.r(rnz) = ri(rnz);
                hr.h(hnz) = hi(hnz);
                dirBoundaryConditions(nodeID) = hr;
            end
            
        end
        
        
        
        function [ Qi, Gi ] = setNeumannBCOnEdge( self, bci )
            %setNeumannBCOnEdge Calculate g and q entries for a Neumann BC
            % This function handles the case when a 'g' and/or 'q' parameter
            % has been included in a pdeBoundaryConditions entity.
            
            
            import pde.internal.elements.*;
            
            N2 = self.N*self.N;
            neqn = self.N*self.numPoints;
            edgeID = bci.appRegionID;
            
            isFuncQ = isa(bci.q, 'function_handle');
            isFuncG = isa(bci.g, 'function_handle');
            anyFunc = isFuncQ || isFuncG;
            anyQ = isFuncQ || any(bci.q(:));
            anyG = isFuncG || any(bci.g(:));
            
            if(~(anyQ || anyG))
                Qi = sparse(neqn,neqn);
                Gi = sparse(neqn,1);
                return;
            end
            [numEdgeNodes, numEdges] = size(bci.elemEdges);
            % preallocate some space for sparse triplets
            preAllocSize = 3*numEdges*self.N;
            if(anyQ)
                qTriplets = zeros(3, preAllocSize);
                numQTriplets = 0;
            end
            if(anyG)
                gTriplets = zeros(2, preAllocSize);
                numGTriplets = 0;
            end
            
            if(numEdgeNodes == 3)
                sfCalc = @shapeLine3;
                dsfCalc = @dShapeLine3;
                numGaussPoints = 3;
            elseif (numEdgeNodes == 2)
                sfCalc = @shapeLine2;
                dsfCalc = @dShapeLine2;
                numGaussPoints = 2;
            else
                error(message('pde:pde2DBCImpl:illegalEdgeDefn', numEdgeNodes));
            end
            
            % compute shape functions and derivatives at gauss points
            intRule=GaussianIntegrationRule(GaussianIntegrationRule.Curve, numGaussPoints);
            intWts = intRule.wts;
            intPoints = intRule.points;
            dsPts = zeros(numEdgeNodes, numGaussPoints);
            for n=1:numGaussPoints
                dsPts(:,n) = dsfCalc(intPoints(:,n));
            end
            sPts = sfCalc(intPoints);
            
            if(anyG)
                Ge = zeros(numEdgeNodes, self.N);
            end
            
            if(anyQ)
                Qe = zeros(numEdgeNodes, numEdgeNodes, N2);
            end
            
            if(~isempty(bci.q) && ~isFuncQ)
                qi = bci.q(:);
                if(any(qi) && length(qi) ~= N2)
                    error(message('pde:pde2DBCImpl:invalidLenQ', edgeID, N2));
                elseif(~any(qi) && length(qi) ~= N2 && length(qi) ~= 1)
                    error(message('pde:pde2DBCImpl:invalidLenQ', edgeID, N2));
                end
            end
            if(~isempty(bci.g) && ~isFuncG)
                gi = bci.g(:);
                if(any(gi) && length(gi) ~= self.N)
                    error(message('pde:pde2DBCImpl:invalidLenG', numEdgeNodes, self.N));
                elseif(~any(gi) && length(gi) ~= self.N && length(gi) ~= 1)
                    error(message('pde:pde2DBCImpl:invalidLenG', numEdgeNodes, self.N));
                end
            end
            
            %  calculate edge normals
            xyAllEdgeNodes = self.points(:, bci.elemEdges(:));
            xyAllEdgeNodes = reshape(xyAllEdgeNodes, 2, numEdgeNodes, numEdges);
            edgeNorm = zeros(numGaussPoints, numEdges);
            elemLength = zeros(numEdges,1);
            if(anyFunc)
                uPoints = zeros(self.N, numGaussPoints, numEdges);
                edgeNormVec = zeros(2, numGaussPoints, numEdges);
            end
            for i=1:numEdges
                xy = xyAllEdgeNodes(:,:,i);
                elemLength(i) = getElemLength(xy(:,1), xy(:,end));
                if(anyFunc)
                    edgeNodes = bci.elemEdges(:,i);
                    uEdgeNodes = self.uN(:, edgeNodes);
                    uPoints(:,:,i) = uEdgeNodes*sPts;
                end
                for j=1:numGaussPoints
                    dxyDr = xy*dsPts(:,j);
                    v2 = cross2([dxyDr;0],[0,0,1]);
                    edgeNorm(j,i) = norm(v2);
                    v2 = v2(1:2);
                    if(anyFunc)
                        edgeNormVec(:,j,i) = v2/edgeNorm(j,i);
                    end
                end
            end
            if(isFuncQ)
                edgeQ = self.callNeumannFuncOnEdge(bci,xyAllEdgeNodes, sPts, bci.q, ...,
                    edgeNormVec, uPoints, 2);
            end
            if(isFuncG)
                edgeG = self.callNeumannFuncOnEdge(bci,xyAllEdgeNodes, sPts, bci.g, ...
                    edgeNormVec, uPoints, 1);
            end
            
            fCol = 1;
            % Iterate over all element edges in this BC
            for i=1:numEdges
                coni = bci.elemEdges(:,i);
                % Iterate over the integration points
                for n=1:numGaussPoints
                    shn = sPts(:,n);
                    detJ = elemLength(i);
                    wt = intWts(n);
                    detWt = detJ*wt;
                    if(anyG)
                        %bci.g
                        if(isFuncG)
                            gi = edgeG(:,fCol);
                        end
                        detShp = detWt*shn;
                        for j=1:self.N
                            Ge(:,j) = Ge(:,j) + gi(j)*detShp;
                        end
                    end
                    if(anyQ)
                        if(isFuncQ)
                            qi = edgeQ(:,fCol);
                        end
                        detShpShp = detWt*(shn*shn');
                        for j=1:N2
                            Qe(:,:,j) = Qe(:,:,j) +qi(j)*detShpShp;
                        end
                    end
                    fCol = fCol + 1;
                end
                if(anyG)
                    % copy the edge g-vector to triplets
                    for j=1:self.N
                        offset = (j-1)*self.numPoints;
                        gj = Ge(:,j);
                        nzInd = gj ~= 0;
                        nnz = sum(nzInd);
                        gTriplets(:, numGTriplets+1:numGTriplets+nnz) = [coni(nzInd)'+offset; gj(nzInd)'];
                        numGTriplets = numGTriplets + nnz;
                    end
                    % set Ge back to zero for the next edge
                    Ge(:) = 0;
                end
                if(anyQ)
                    % copy the edge q-matrix to triplets
                    jj = 1;
                    for j1=1:self.N
                        offset1 = (j1-1)*self.numPoints;
                        for j2=1:self.N
                            offset2 = (j2-1)*self.numPoints;
                            qjk = reshape(Qe(:, :, jj)', [], 1);
                            nzInd = qjk ~= 0;
                            nnz = sum(nzInd);
                            if(nnz)
                                gRows = repmat(coni' + offset2, numEdgeNodes, 1); gRows = gRows(:);
                                gCols = repmat(coni' + offset1, 1, numEdgeNodes);
                                qTriplets(:,numQTriplets+1:numQTriplets + nnz) = [gRows(nzInd)'; gCols(nzInd); qjk(nzInd)'];
                                numQTriplets = numQTriplets + nnz;
                            end
                            jj = jj + 1;
                        end
                    end
                    % set Qe back to zero for the next edge
                    Qe(:) = 0;
                end
            end
            
            if(anyG && numGTriplets)
                Gi = sparse(gTriplets(1, 1:numGTriplets), 1, gTriplets(2, 1:numGTriplets), ...
                    neqn,1);
            else
                Gi = sparse(neqn,1);
            end
            if(anyQ && numQTriplets)
                Qi = sparse(qTriplets(1, 1:numQTriplets), qTriplets(2, 1:numQTriplets), ...
                    qTriplets(3, 1:numQTriplets), neqn, neqn);
            else
                Qi = sparse(neqn,neqn);
            end
            
        end
        
        
        function setDirichletBCOnEdge( self, bci, dirBoundaryConditions )
            %setDirichletBCOnEdge Calculate r and h entries for a Dirichlet BC
            % This function handles the case when a 'r' and/or 'h' parameter has been
            % included in a pdeBoundaryConditions entity.
            
            import pde.internal.*;
            
            edgeID = bci.appRegionID;
            numNodes = length(bci.nodes);
            
            % Handle the first term (r vector)
            isFuncR = isa(bci.r, 'function_handle');
            isFuncH = isa(bci.h, 'function_handle');
            
            nodesCell = num2cell(bci.nodes);
            doesNodeHaveBC = dirBoundaryConditions.isKey(nodesCell);
            
            
            % Following flags are used to check if any node on this FE edge
            % or any DoF at a node still needs BC assignment. If not, do
            % nothing and return to the caller.
            
            nodesNeedBCAssignment = false;
            dofsNeedBCAssignment = false;
            if (~all(doesNodeHaveBC))
                nodesNeedBCAssignment = true;
            end
            
            if (nodesNeedBCAssignment == false)
                % check for BC at all DoFs only when all nodes already have some BC
                
                for n = cell2mat(nodesCell(doesNodeHaveBC))'
                    hr = dirBoundaryConditions(n);
                    previousBCAssignment = hr.r ~=0;
                    for j = 1:numel(hr.r)
                        if (~previousBCAssignment(j)) % If there was no previous assignment
                            dofsNeedBCAssignment = true;
                            continue
                        end
                    end
                end
            end
            
            if ~(nodesNeedBCAssignment || dofsNeedBCAssignment)
                return
            end
            
            if(isFuncR || isFuncH)
                edgeNormals = self.calcEdgeNormalsAtNodes(bci.elemEdges, bci.nodes);
            end
            if(~isempty(bci.term1))
                if isFuncR
                    rMat = self.callDirichletFuncOnEdge(bci, bci.nodes, bci.r, edgeNormals, 1);
                else
                    ri = bci.term1(:);
                    if(length(ri) ~= self.N)
                        error(message('pde:pde2DBCImpl:invalidLenR', edgeID, self.N));
                    end
                    rnz = ri~=0;
                end
            else
                ri = zeros(self.N,1);
                rnz = ri~=0;
            end
            
            % Handle the second term (h matrix)
            if(isempty(bci.term2))
                hi = eye(self.N);
                hnz = hi~=0;
            else
                if isFuncH
                    hMat = self.callDirichletFuncOnEdge(bci, bci.nodes, bci.h, edgeNormals, 2);
                    hMat = reshape(hMat, self.N, self.N, []);
                else
                    hi = bci.term2(:);
                    if(length(hi) ~= self.N^2)
                        error(message('pde:pde2DBCImpl:invalidLenH', edgeID, self.N^2));
                    end
                    hnz = hi~=0;
                end
            end
            
            hrZ = HRNode(self.N);
            
            
            for i=1:numNodes
                nodeID = bci.nodes(i);
                if(doesNodeHaveBC(i))
                    hr = dirBoundaryConditions(nodeID);
                else
                    hr = hrZ;
                end
                if(isFuncH)
                    hi = hMat(:,:,i);
                    hnz = hi~=0;
                end
                if(isFuncR)
                    ri = rMat(:,i);
                    rnz = ri~=0;
                end
                hr.r(rnz) = ri(rnz);
                hr.h(hnz) = hi(hnz);
                dirBoundaryConditions(nodeID) = hr;
            end
        end
        
        
        
        function rh = callDirichletFuncOnEdge( self, bc, nodes, func, edgeNormals, t12 )
            %callDirichletFuncOnEdge Call a user-defined function for a Dirichlet BC
            % This function handles the case when a function handle has been defined
            % for the 'h' or 'r' parameter in a pdeBoundaryConditions entity.
            
            
            % If t12 is 1, this function is being called for an r-vector, otherwise
            % and h-matrix
            if(t12 == 1)
                numRows = self.N;
            else
                numRows = self.N^2;
            end
            n = length(nodes);
            p = self.points;
            if(self.gotTime)
                state.time = self.time;
            end
            if(bc.isVectorized)
                appRegion = self.applicationRegion(p(:,nodes), edgeNormals);
                state.u = self.uN(:,nodes);
                if nargin(func) == 3
                    rh = func(self.myPde, appRegion, state);
                else
                    rh = func(appRegion, state);
                end
                if(t12 == 1)
                    self.checkFuncEvalVecN('Dirichlet', func, rh, n);
                else
                    self.checkFuncEvalVecNxN('Dirichlet', func, rh, n);
                end
                rh = reshape(rh, [], n);
            else
                xyPts = p(:,nodes);
                rh = zeros(numRows, n);
                for i=1:n
                    ni = nodes(i);
                    appRegion = self.applicationRegion(xyPts(:,i), edgeNormals(:,i));
                    state.u = self.uN(:,ni);
                    if nargin(func) == 3
                        bci = func(self.myPde, appRegion, state);
                        bci = bci(:);
                    else
                        bci = func(appRegion, state);
                        bci = bci(:);
                    end
                    if(length(bci) ~= numRows)
                        error(message('pde:pde2DBCImpl:invalidNumEntriesDirFunc', func2str(func), numRows));
                    end
                    rh(:,i) = bci(:);
                end
            end
            
        end
        
        
        function bc = callNeumannFuncOnEdge(self, bc, xy, sPts, func, edgeNormVec, ...
                uPoints, t12)
            %callNeumannFuncOnEdge Call a user-defined function for Neumann BC
            % This function handles the case when a function handle has been defined
            % for the 'g' or 'h' parameter in a pdeBoundaryConditions entity.
            
            
            if(t12 == 1)
                numRows = self.N;
            else
                numRows = self.N^2;
            end
            numEdges = size(xy, 3);
            numIntPts = size(sPts, 2);
            xyPts = zeros(2, numIntPts, numEdges);
            if(self.gotTime)
                state.time = self.time;
            end
            for i=1:numEdges
                xyPts(:,:,i) = xy(:,:,i)*sPts;
            end
            xyPts = reshape(xyPts, 2, []);
            uPoints = reshape(uPoints, self.N, []);
            edgeNormVec = reshape(edgeNormVec, 2, []);
            n = size(xyPts,2);
            if(bc.isVectorized && n>1)
                % If the user's function can handle vector arguments and
                % we need to evaluate at more than one point, call the function
                % once.
                appRegion = self.applicationRegion(xyPts, edgeNormVec);
                state.u = uPoints;
                if nargin(func) == 3
                    bc = func(self.myPde, appRegion, state);
                else
                    bc = func(appRegion, state);
                end
                if(t12 == 1)
                    self.checkFuncEvalVecN('Neumann', func, bc, n);
                else
                    self.checkFuncEvalVecNxN('Neumann', func, bc, n);
                end
                bc = reshape(bc, [], n);
            else
                bc = zeros(numRows, n);
                for i=1:n
                    appRegion = self.applicationRegion(xyPts(:,i), edgeNormVec(:,i));
                    state.u = uPoints(:,i);
                    if nargin(func) == 3
                        bci = func(self.myPde, appRegion, state);
                    else
                        bci = func(appRegion, state);
                    end
                    bci = bci(:);
                    if(length(bci) ~= numRows)
                        error(message('pde:pde2DBCImpl:invalidNumEntriesDirFunc', func2str(func), numRows));
                    end
                    bc(:,i) = bci;
                end
            end
            
        end
        
        
        
        function rMat = callValueFuncOnEdge(self, bc, nodes, func, edgeNormals)
            %callValueFuncOnEdge Call user-defined function for the value vector
            % The vector of values in a BC of type Value can be defined in a user-written
            % function.
            
            n = length(nodes);
            p = self.points;
            if(self.gotTime)
                state.time = self.time;
            end
            if(bc.isVectorized)
                appRegion = self.applicationRegion(p(:,nodes), edgeNormals);
                state.u = self.uN(:,nodes);
                if nargin(func) == 3
                    rMat = func(self.myPde, appRegion, state);
                else
                    rMat = func(appRegion, state);
                end
                numRet = size(rMat, ndims(rMat));
                if(numRet ~= n)
                    error(message('pde:pde2DBCImpl:invalidNumDirFunc', func2str(func), n, numRet));
                end
            else
                % If the user's function is not vectorized, call it for each edge.
                xyPts = p(:,nodes);
                for i=1:n
                    ni = nodes(i);
                    appRegion = self.applicationRegion(xyPts(:,i), edgeNormals(:,i));
                    state.u = self.uN(:,ni);
                    if nargin(func) == 3
                        bci = func(self.myPde, appRegion, state);
                    else
                        bci = func(appRegion, state);
                    end
                    if(i == 1)
                        rMat = zeros(length(bci), n);
                    end
                    rMat(:,i) = bci(:);
                end
            end
            
            
        end
        
        
        
        
        function edgeNormVecAtNodes = calcEdgeNormalsAtNodes(self, edges, nodes)
            %calcEdgeNormalsAtNodes Calculate normal to a edge at each of its nodes
            % Given a set of edges and edge nodes, this function calculates
            % the edge normal at each node.
            
            import pde.internal.elements.*;
            
            [numEdgeNodes, numEdges] = size(edges);
            
            % Get the shape function derivatives wrt the natural coordinates
            % for the line shape that represents a edge of a triangle.
            if(numEdgeNodes == 3)
                dsfCalc = @dShapeLine3;
                rs = [-1,0,1]'; % natural coordinates at the 3 nodes
                numNodesToEval = 3;
            elseif (numEdgeNodes == 2)
                dsfCalc = @dShapeLine2;
                rs = [-1,1]'; % natural coordinates at the 2 nodes
                numNodesToEval = 2;
            else
                error(message('pde:pde2DBCImpl:illegalEdgeDefn', numEdgeNodes));
            end
            % Collect the shape function derivatives for all edge nodes
            dsNodes = zeros(numEdgeNodes, numEdgeNodes);
            for n=1:numEdgeNodes
                dsNodes(:,n) = dsfCalc(rs(n));
            end
            
            numNodes = length(nodes);
            numEdgesAtNode = zeros(numNodes, 1);
            edgeNormVecAtNodes = zeros(2, numNodes);
            
            nodeIndexVec = zeros(self.numPoints, 1);
            nodeIndexVec(nodes) = 1:numNodes;
            
            xyAllEdgeNodes = self.points(:, edges(:));
            xyAllEdgeNodes = reshape(xyAllEdgeNodes, 2, numEdgeNodes, numEdges);
            edgeNormVec = zeros(2, numEdgeNodes);
            for i=1:numEdges
                xy = xyAllEdgeNodes(:,:,i);
                for j=1:numNodesToEval
                    % Compute the vectors dX/dr
                    dxyDr = xy*dsNodes(:,j);
                    % Their cross product is the normal at the point
                    v2 = cross2([dxyDr;0],[0,0,1]);
                    v2 = v2(1:2);
                    edgeNormVec(:,j) = v2/norm(v2);
                end
                if(numNodesToEval==1)
                    edgeNormVec = repmat(edgeNormVec(:,1), 1, numEdgeNodes);
                end
                for j=1:numEdgeNodes
                    nj = nodeIndexVec(edges(j,i));
                    numEdgesAtNode(nj) = numEdgesAtNode(nj) + 1;
                    edgeNormVecAtNodes(:,nj) = edgeNormVecAtNodes(:,nj) + edgeNormVec(:,j);
                end
            end
            
            % The normal at a node is the average of the normals of all the edges
            % connected to that node.
            for i=1:numNodes
                edgeNormVecAtNodes(:,i) = edgeNormVecAtNodes(:,i)/numEdgesAtNode(i);
            end
            
        end
        
        
        function de = extractDiscreteEdgeQuadratic(~,e,eid)
            % Given the e from a p, e, t mesh, extract the discrete edge
            % corresponding to the geometry edge whos ID is eid.
            
            % The e matrix in 2-D treats two segments of a quadratic
            % element as two separate edges. The following code stiches
            % adjacent edges to obtain a quadratic edge. Domain to the left
            % or right of edge info, available in e matrix, is used to
            % preserve orientation of quadratic edges.
            
            idx = (e(5,:) == eid);
            subDomainToLeft = e(6, idx);
            
            facets = e([1, 2], idx);
            numfacets = size(facets,2);
            oddidx = 1:2:numfacets;
            evenidx = 2:2:numfacets;
            de = [facets(:,oddidx);facets(2,evenidx)];
            sizeDe = size(de);
            if (~all(subDomainToLeft))
                de = flipud(de(:));
            end
            de = reshape(de,sizeDe);
        end
        
        function de = extractDiscreteEdgeLinear(~,e,eid)
            % Given the e from a p, e, t mesh, extract the discrete edge
            % corresponding to the geometry edge whos ID is eid.
            
            %Unlike quadratic elements, no further manipulation is needed
            %in linear element case.
            idx = (e(5,:) == eid);
            de = e([1, 2], idx);
        end
        
    end % methods
    
    
    methods(Static, Access=private)
        function nodes=getNodesForEdges(edges)
            nodes = unique(edges(:));
        end
        
        function appRegion = applicationRegion(locations, edgeNormals )
            %APPLICATIONREGION Construct application region struct for user function
            appRegion.x = locations(1,:);
            appRegion.y = locations(2,:);
            appRegion.nx = edgeNormals(1,:);
            appRegion.ny = edgeNormals(2,:);
        end
    end % methods
    
    
end


function hN = getH(N, dofList, edgeID)
% get indices into  an N*N h-vector for the list of dofs
lenDofs = length(dofList);
hN = zeros(N);
for k=1:lenDofs
    kk = dofList(k);
    if(kk < 1 || kk > N || ~isreal(kk) || mod(kk,1) ~= 0)
        kkStr = num2str(kk, 16);
        error(message('pde:pde2DBCImpl:eqnNumInvalidValue', kkStr, edgeID, N));
    end
    hN(kk,kk) = 1;
end
end

function checkTermLengths(lenT1, lenT2, edgeID)
if(lenT1~=lenT2)
    eId = 'pde:pde2DBCImpl:valUnequalLengths';
    err = message(eId, lenT1, lenT1, lenT2, edgeID);
    throwAsCaller(MException(eId, err.getString()));
end
end


% This highly specialized version of the cross function with no error
% checking is used to improve performance.
function c=cross2(a,b)
c=[a(2)*b(3)-a(3)*b(2);
    a(3)*b(1)-a(1)*b(3);
    a(1)*b(2)-a(2)*b(1)];
end

function elmLength=getElemLength(p1,p2)
% compute length of a line segment given two end point coordinates.
elmLength = sqrt((p2(1)-p1(1))^2 + (p2(2)-p1(2))^2);

end