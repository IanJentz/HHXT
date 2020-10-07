function [results] = calcGradients(results,model)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Function wide variables

    n1 = model.PDESystemSize;
    matlDefs = model.MaterialProperties.MaterialAssignments;
    matlDefs2 = copy(model.MaterialProperties.MaterialAssignments);
    matlMap = model.MaterialMap;
    mesh = model.Mesh;
    [~,nNodes] = size(mesh.Nodes);
    order = model.Mesh.GeometricOrder;
    
    analysistype = model.AnalysisType;
        
    if results.IsTwoD
        nRegions = model.Geometry.NumFaces;
        nBoundaries = model.Geometry.NumEdges;
        regionType = 'Face';
        boundaryType = 'Edge';
        switch order
            case 'linear'
                element_type = 'MatLab3NodeTriangle';
            case 'quadratic'
                element_type = 'MatLab6NodeTriangle';
        end
    else % 3D
        nRegions = model.Geometry.NumCells;
        nBoundaries = model.Geometry.NumFaces;
        regionType = 'Cell';
        boundaryType = 'Face';
        switch order
            case 'linear'
                element_type = 'MatLab4NodeTet';
            case 'quadratic'
                element_type = 'MatLab10NodeTet';
        end
    end
    
    times = results.SolutionTimes;
    nT = length(times);
    time_ind = 1:nT;
    state = result2state(results,time_ind);
    
%% Calculate Gradients    
    rst = nodeRST('all',element_type);
    [~,npere] = size(rst);

    % initialize derivatives
    dudx = zeros(nNodes,n1,nT); dudy = dudx;
    if ~results.IsTwoD; dudz = dudx; end
    divby = zeros(nNodes,n1,nT); % keep track of how many elements were queried for gradient calculation
    
    % streams
    l1 = (n1-1)/2;
    streams = 1:l1;     
    n2l = repmat(1:l1,2,1)+1; n2l = [1;n2l(:)];
    
    for cellID = 1:nRegions
    e = findElements(mesh,'region',regionType,cellID);
    tm = mesh.Elements(:,e);
    [~,ne] = size(tm);

    for k = 1:ne
        pe = tm(:,k); % element nodes
        X = mesh.Nodes(:,pe)'; % node positions
        unodal = state.u(:,pe,:);
        
        for i = 1:npere
            %[N,B,J,detJ] = shapeFunction(rst,X,element_type)
            [~,B,~,detJ] = shapeFunction(rst(:,i),X,element_type); 

            for n = 1:n1
                % determine stream and PDE position n
                l = n2l(n);
                
                
                if matlMap(l,cellID) ~= 0
                for ti = 1:nT
                    DELu = detJ*B*unodal(n,:,ti)';

                    dudx(pe(i),n,ti) = dudx(pe(i),n) + DELu(1); 
                    dudy(pe(i),n,ti) = dudy(pe(i),n) + DELu(2);
                    if ~results.IsTwoD
                        dudz(pe(i),n,ti) = dudz(pe(i),n) + DELu(3);
                    end
                    divby(pe(i),n,ti) = divby(pe(i),n) + detJ; 
                end
                end
            end               

        end % loop over element nodes
    end     

    end % loop over regions

    % divide by the number of elements used to calculate the gradient
    dudx = dudx./divby;
    dudy = dudy./divby;
    if ~results.IsTwoD; dudz = dudz./divby; end

    % check for NaN results       
    dudx(isnan(dudx)) = 0;
    dudy(isnan(dudy)) = 0;
    if ~results.IsTwoD; dudz(isnan(dudz)) = 0; end

    for n = 1:n1
        results.XGradients(:,n,:) = dudx(:,n,:);
        results.YGradients(:,n,:) = dudy(:,n,:);
        if ~results.IsTwoD; results.ZGradients(:,n,:) = dudz(:,n,:); end

    end
        
    
%% Internal Functions
function rst = nodeRST(node_index,element_type)

    switch element_type
    
    case 'MatLab10NodeTet'
    %  Matlab has node order like this:
    %                    4
    %                    o
    %                   /|\_
    %                  / |  \_
    %                 /  |    \_
    %                /   o 8    \o 10
    %               /    |        \_
    %            9 o     |          \_
    %             /      | 1     7    \_
    %            /       o-------o-------o 3
    %           /   5 __/         _____/
    %          /   _o/      _____/ 
    %         / __/   _____/o   
    %        /_/_____/       6
    %     2 o__/	 

    % return rst of each point
    switch node_index
        case 1
            rst = [0;0;0];
        case 2
            rst = [1;0;0];
        case 3
            rst = [0;1;0];
        case 4
            rst = [0;0;1];
        case 5
            rst = [0.5;0;0;];
        case 6
            rst = [0.5;0.5;0];
        case 7
            rst = [0;0.5;0];
        case 8
            rst = [0;0;0.5];
        case 9
            rst = [0.5;0;0.5];
        case 10
            rst = [0;0.5;0.5];
        case 'all'
            rst = [ 0 , 1 , 0 , 0 ,0.5,0.5, 0 , 0 ,0.5 ;...
                    0 , 0 , 1 , 0 , 0 ,0.5,0.5, 0 , 0 ;...
                    0 , 0 , 0 , 1 , 0 , 0 , 0 ,0.5, 0.5 ];
    end
        
    case {'MatLab4NodeTet','SOLID70'}
    %  Matlab has node order like this:
    %                    4
    %                    o
    %                   /|\_
    %                  / |  \_
    %                 /  |    \_
    %                /   |      \_ 
    %               /    |        \_
    %              /     |          \_
    %             /      | 1          \_
    %            /       o---------------o 3
    %           /     __/         _____/
    %          /   __/      _____/ 
    %         / __/   _____/   
    %        /_/_____/       
    %     2 o__/	 

    % return rst of each point
    switch node_index
        case 1
            rst = [0;0;0];
        case 2
            rst = [1;0;0];
        case 3
            rst = [0;1;0];
        case 4
            rst = [0;0;1];
        case 'all'
            rst = [ 0 , 1 , 0 , 0  ;...
                    0 , 0 , 1 , 0  ;...
                    0 , 0 , 0 , 1  ];
    end

    case 'MatLab6NodeTriangle'
    %  Matlab has node order like this:
    %      
    %   3 o 
    %     |\ 
    %     | \
    %     |  \
    %     |   \
    %   6 o    o 5
    %     |     \
    %     |      \
    %     |       \
    %     |        \
    %   1 o----o----o 2
    %          4

    % return rst of each point
    switch node_index
        case 1
            rst = [0;0];
        case 2
            rst = [1;0];
        case 3
            rst = [0;1];
        case 4
            rst = [0.5;0];
        case 5
            rst = [0.5;0.5];
        case 6
            rst = [0;0.5];
        case 'all'
            rst = [ 0 , 1 , 0 ,0.5,0.5, 0  ;...
                    0 , 0 , 1 , 0 ,0.5,0.5 ];
    end
    
    case 'MatLab3NodeTriangle'
    %  Matlab has node order like this:
    %      
    %   3 o 
    %     |\ 
    %     | \
    %     |  \
    %     |   \
    %     |    \
    %     |     \
    %     |      \
    %     |       \
    %     |        \
    %   1 o---------o 2
    %          
    
    % return rst of each point
    switch node_index
        case 1
            rst = [0;0];
        case 2
            rst = [1;0];
        case 3
            rst = [0;1];
        case 'all'
            rst = [ 0 , 1 , 0 ;...
                    0 , 0 , 1 ];
    end
    

    end
    
end


function [N,B,J,detJ] = shapeFunction(rst,X,element_type)
%SHAPEFUNCTION returns B and J matrixes of an element with nodal locations
% X at dimensionless element point rst = [r, s, t]
%   [N,B,J,detJ] = shapeFunction(rst,X,element_type) returns matrices of
%   element of element_type with nodal locations X at dimensionless element
%   point rst = [r,s,t].
%
%   element_type: 'MatLab10NodeTet','SOLID187' - 10 node tetrahedron
%                 'MatLab4NodeTet','SOLID70' - 4 node tetrahedron
%                 'MatLab6NodeTriangle' - 6 node triangle
%                 'MatLab3NodeTriangle' - 4 node triangle
    justN = false;
    if nargout == 1
        justN = true;
    end
    % calculate the shape function partial derivatives based on the element
    % type
    switch element_type
        
    case 'SOLID187'
    %  Plesha/Witt Fig 7.1-4 has node order like this:
    %                    4
    %                    o
    %                   /|\_
    %                  / |  \_
    %                 /  |    \_
    %                /   o 7    \o 9
    %               /    |        \_
    %           10 o     |          \_
    %             /      | 1     6    \_
    %            /       o-------o-------o 3
    %           /   5 __/         _____/
    %          /   _o/      _____/ 
    %         / __/   _____/o   
    %        /_/_____/       8
    %     2 o__/	 

    % take in evaluation points
    r = rst(1); s = rst(2); t = rst(3);
    % calculate the shape functions
    N = [ (1-r-2-t)*(1-2*r-2*s-2*t),...
          r*(2*r-1),s*(2*s-1),t*(2*t-1),...
          4*r*(1-r-s-t),4*s*(1-r-s-t),4*t*(1-r-s-t),...
          4*rs*s,4*s*t,4*t*r ];
    % calculate the shape function partial derivatives for a 10 node Tet
    rst1 = (r+s+t)-0.75;    % used for values in column 1 of p_N
    rst2 = (1-r-s-t);       % used for vlaues in coulmns 5-7 of p_N
    if justN == false
    % as per element in Plesh/Witt Fig 7.1-4
    p_N = 4*...
    [ rst1, r-.25,     0,     0, rst2-r,     -s,     -t, s, 0, t ;...
      rst1,     0, s-.25,     0,     -r, rst2-s,     -t, r, t, 0 ;...
      rst1,     0,     0, t-.25,     -r,     -s, rst2-t, 0, s, r ];
    % use first 10 nodes
    X = X(1:10,:);
    end
    
    case 'MatLab10NodeTet'
    %  Matlab has node order like this:
    %                    4
    %                    o
    %                   /|\_
    %                  / |  \_
    %                 /  |    \_
    %                /   o 8    \o 10
    %               /    |        \_
    %            9 o     |          \_
    %             /      | 1     7    \_
    %            /       o-------o-------o 3
    %           /   5 __/         _____/
    %          /   _o/      _____/ 
    %         / __/   _____/o   
    %        /_/_____/       6
    %     2 o__/	 

    % take in evaluation points
    r = rst(1); s = rst(2); t = rst(3);
    % calculate the shape functions
    N = [ (1-r-2-t)*(1-2*r-2*s-2*t),...
          r*(2*r-1),s*(2*s-1),t*(2*t-1),...
          4*r*(1-r-s-t),...
          4*t*(1-r-s-t),4*rs*s,...
          4*s*(1-r-s-t),4*t*r,4*s*t ];
    % calculate the shape function partial derivatives for a 10 node Tet
    rst1 = (r+s+t)-0.75;    % used for values in column 1 of p_N
    rst2 = (1-r-s-t);       % used for vlaues in coulmns 5-7 of p_N
    if justN == false
    % as per MatLab's 'quadratic' 3D (10 node tet) element
    p_N = 4*...
    [ rst1, r-.25,     0,     0, rst2-r,     -t, s,     -s, t, 0 ;...
      rst1, 0    , s-.25,     0,     -r,     -t, r, rst2-s, 0, t ;...
      rst1, 0    ,     0, t-.25,     -r, rst2-t, 0,     -s, r, s ];
    % use first 10 nodes
    X = X(1:10,:);
    end
    
    case {'MatLab4NodeTet','SOLID70'}
    %  Matlab has node order like this:
    %                    4
    %                    o
    %                   /|\_
    %                  / |  \_
    %                 /  |    \_
    %                /   |      \_ 
    %               /    |        \_
    %              /     |          \_
    %             /      | 1          \_
    %            /       o---------------o 3
    %           /     __/         _____/
    %          /   __/      _____/ 
    %         / __/   _____/   
    %        /_/_____/       
    %     2 o__/	 

    % take in evaluation points
    r = rst(1); s = rst(2); t = rst(3);
    % calculate the shape functions
    N = [ (1-r-s-t),r,s,t ];
    if justN == false
    % 4 node tet element
    p_N = 1*...
    [ -1, 1, 0, 0 ;...
      -1, 0, 1, 0 ;...
      -1, 0, 0, 1 ];
    % use first 4 nodes
    X = X(1:4,:);
    end
    
  case 'MatLab6NodeTriangle'
    %  Matlab has node order like this:
    %      
    %   3 o 
    %     |\ 
    %     | \
    %     |  \
    %     |   \
    %   6 o    o 5
    %     |     \
    %     |      \
    %     |       \
    %     |        \
    %   1 o----o----o 2
    %          4

    % take in evaluation points
    r = rst(1); s = rst(2);
    % calculate the shape functions
    N = [ (1-r-s)*(1-2*r-2*s),...
           r*(2*r-1),s*(2*s-1),...
           4*r*(1-r-s),4*r*s,4*s*(1-r-s)];
    if justN == false
    % 6 node triangle element
    r2s = 2*r+s;
    rs2 = r+2*s;
    p_N = 4*...
    [ r2s-0.75, r-0.25,      0, 1-r2s, s, -s    ;...
      rs2-0.75,      0, s-0.25,    -r, r, 1-rs2 ];
    % use first 6 nodes
    X = X(1:6,:);
    end
    
    case 'MatLab3NodeTriangle'
    %  Matlab has node order like this:
    %      
    %   3 o 
    %     |\ 
    %     | \
    %     |  \
    %     |   \
    %     |    \
    %     |     \
    %     |      \
    %     |       \
    %     |        \
    %   1 o---------o 2
    %          
    
    % take in evaluation points
    r = rst(1); s = rst(2);
    % calculate the shape functions
    N = [ (1-r-s),r,s ];
    if justN == false
    % 3 node triangle element
    p_N = 1*...
    [ -1, 1, 0 ;...
      -1, 0, 1 ];
    % use first 3 nodes
    X = X(1:3,:);
    end
    
    end
    
    if justN == false
    % calculate the Jacobian and the B matrix
    J = p_N*X;
    detJ = det(J);
    B = J\p_N;
    end

end    
    
end

