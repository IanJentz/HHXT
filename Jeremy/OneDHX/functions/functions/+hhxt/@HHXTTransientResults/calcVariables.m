function [results] = calcVariables(results,model)
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
    
 %% establish state (single time)
    
times = results.SolutionTimes;
nT = length(times);
time_ind = 1:nT;
state = result2state(results,time_ind);
    
%% Parse Coefficient Matrix Options
varargin = model.MatrixOptions;

resistancemode = 'resistance';
Rfc_mult = 1;

for i = 1:2:length(varargin)
    switch varargin{i}
        case {'resistancemode','ResistanceMode'}
            switch varargin{i+1}
                case {'fixed','resistance'}
                    resistancemode = 'resistance';
                case {'heattransfer','colburn','Colburn'}
                    resistancemode = 'colburn';
                case {'mixed','Mixed'}
                    resistancemode = 'mixed';
                otherwise
                    error(['unrecognised input: ',varargin{i},',',varargin{i+1}])
            end 
        case {'ResistanceMultiplier'}
            Rfc_mult = varargin{i+1};
    end

end
    
%% Variables;

    % streams
    l1 = (n1-1)/2;
    streams = 1:l1;

    % initialize variables
    uDm = zeros(nNodes,l1,nT); uDx = uDm; uDy = uDm; if ~results.IsTwoD; uDz = uDm; end
    uCm = zeros(nNodes,l1,nT); uCx = uCm; uCy = uCm; if ~results.IsTwoD; uCz = uCm; end
    Rem = zeros(nNodes,l1,nT); Rex = Rem; Rey = Rem; if ~results.IsTwoD; Rez = Rem; end
    sH = zeros(nNodes,l1,nT); sR = sH;
    
    divby = zeros(nNodes,l1,nT); % keep track of how many elements were queried for gradient calculation
    
    for ti = 1:nT
        
        time = state.time(ti);
        
        uFull = state.u(:,:,ti)';
        uFull = uFull(:);
        % update permeability
        tempresult = hhxt.HHXTResults(model,uFull,0);
        model2 = updatePermeability(model,tempresult);
        model2.HHXTSolverOptions.curtime = time;
        matlDefs = model2.MaterialProperties.MaterialAssignments;
    
    for cellID = 1:nRegions
        j = findNodes(mesh,'region',regionType,cellID);
        pos = mesh.Nodes(:,j);
        
        for str = 1:l1
            % determine stream and PDE position n
            l = streams(str)+1;
            n = 2*l-2;
            
            if matlMap(l,cellID) ~= 0
            
            % Resitance terms applying to u(l)---------------------------------
        
            % determine properties
            phi = matlDefs(matlMap(l,cellID)).PhaseFraction;
            rho = matlDefs(matlMap(l,cellID)).MassDensity;
            mu = matlDefs(matlMap(l,cellID)).Viscosity;
            Dh = matlDefs(matlMap(l,cellID)).HydraulicDiameter;    
            k = matlDefs(matlMap(l,cellID)).ThermalConductivity;
            cp = matlDefs(matlMap(l,cellID)).SpecificHeat;
            if isa(phi,'function_handle') || isa(rho,'function_handle') ...
                   || isa(mu,'function_handle') || isa(Dh,'function_handle')...
                   || isa(k,'function_handle') || isa(cp,'function_handle')
                % determine temp, pressure, and pressure gradient
                temp_fluid = state.u(n+1,j,ti); %20;  
                press = state.u(n,j,ti);
            end
            if isa(phi,'function_handle')
                phi = phi(pos,time);
            end
            if isa(rho,'function_handle')
                rho = rho(pos,time,temp_fluid,press);
            end
            if isa(mu,'function_handle')
                mu = mu(pos,time,temp_fluid,press);
            end
            if isa(Dh,'function_handle')
                Dh = Dh(pos,time);
            end 
            if isa(k,'function_handle')
                k = k(pos,time,temp_fluid,press);
            end
            if isa(cp,'function_handle')
                cp = cp(pos,time,temp_fluid,press);
            end
            
%             % interpolate permiability from nodal result maps
%             k_results = matlDefs(matlMap(l,cellID)).NodalPermeability;
%             k_perm = interpolateNodalVector(k_results,pos,isotropicPermeability);
            
            % determine Darcy velocity and Reynolds number
            % note u = uD/phi
%             uD = -k_perm'.*DELP; % note u = u_D/phi;    
            uD_results = matlDefs(matlMap(l,cellID)).NodalDarcyVelocity;
            uD = uD_results.NodalSolution(j,:);

            if results.IsTwoD
                Re = (rho.*Dh./(mu.*phi))'.*uD(:,[1,1]);
            else
                Re = (rho.*Dh./(mu.*phi))'.*uD(:,[1,1,1]);
            end
            
            uC = uD./phi; % channel fluid velocity

            uDm(j,l-1,ti) = uDm(j,l-1,ti) + sqrt(sum(uD.^2,2));
            uDx(j,l-1,ti) = uDx(j,l-1,ti) + uD(:,1);
            uDy(j,l-1,ti) = uDy(j,l-1,ti) + uD(:,2);
            if ~results.IsTwoD; uDz(j,l-1,ti) = uDz(j,l-1,ti) + uD(:,3); end
            uCm(j,l-1,ti) = uCm(j,l-1,ti) + sqrt(sum(uC.^2,2));
            uCx(j,l-1,ti) = uCx(j,l-1,ti) + uC(:,1);
            uCy(j,l-1,ti) = uCy(j,l-1,ti) + uC(:,2);
            if ~results.IsTwoD; uCz(j,l-1,ti) = uCz(j,l-1,ti) + uC(:,3); end
            Rem(j,l-1,ti) = Rem(j,l-1,ti) + sqrt(sum(Re.^2,2));
            Rex(j,l-1,ti) = Rex(j,l-1,ti) + Re(:,1);
            Rey(j,l-1,ti) = Rey(j,l-1,ti) + Re(:,2);
            if ~results.IsTwoD; Rez(j,l-1,ti) = Rez(j,l-1,ti) + Re(:,3); end
            
            
            uD = sqrt(sum(uD.*uD,2));
            Re = sqrt(sum(Re.*Re,2));
             
            
            % determine Resistance
            switch resistancemode 
                case 'resistance'
                    % determine resistance between stream and wall
                    Rfc = matlDefs(matlMap(l,cellID)).Resistance;
                    if isa(Rfc,'function_handle')
                        Rfc = Rfc(pos,time,Re); % resitance between fluid and core can be position and temperature dependent
                    end
                case 'colburn'
                    if ~isempty(matlDefs(matlMap(l,cellID)).Resistance)
                        % determine resistance between wall and core
                        Rfc = matlDefs(matlMap(l,cellID)).Resistance;
                        if isa(Rfc,'function_handle')
                            Rfc = Rfc(pos,time,Re); % resitance between fluid and core can be position and temperature dependent
                        end
                    else
                        Rfc = 0;
                    end
                    % determine colburn heat transfer coefficient
                    jC = matlDefs(matlMap(l,cellID)).jColburn;
                    if isa(jC,'function_handle')
                        jC = jC(pos,time,Re);
                    end
                    % determine resistance
                    Rfc = ( (4./Dh).*jC.*(rho./phi).*cp.^(1/3).*k.^(2/3).*mu.^(-2/3).*uD ).^(-1) + Rfc;
                case 'mixed'
                    if ~isempty(matlDefs(matlMap(l,cellID)).Resistance)
                        % determine resistance between wall and core
                        Rfc = matlDefs(matlMap(l,cellID)).Resistance;
                        if isa(Rfc,'function_handle')
                            Rfc = Rfc(pos,time,Re); % resitance between fluid and core can be position and temperature dependent
                        end
                    else
                        Rfc = 0;
                    end
                    if ~isempty(matlDefs(matlMap(l,cellID)).jColburn)
                        % determine colburn heat transfer coefficient
                        jC = matlDefs(matlMap(l,cellID)).jColburn;
                        if isa(jC,'function_handle')
                            jC = jC(pos,time,Re);
                        end
                        % determine resistance
                        Rfc = ( (4./Dh).*jC.*(rho./phi).*cp.^(1/3).*k.^(2/3).*mu.^(-2/3).*uD ).^(-1) + Rfc;  
                    end    
            end % switch resistance mode
            
            Rfc = Rfc*Rfc_mult; % volumetric thermal resistance
                        
            Qvol = ( state.u(1,j,ti)-state.u(n+1,j,ti) )'./Rfc; % volumetric heat flowing into stream
            
            sH(j,l-1,ti) = sH(j,l-1,ti) + Qvol(:,1);
            sR(j,l-1,ti) = sR(j,l-1,ti) + Rfc(:,1);
            
            divby(j,l-1,ti) = divby(j,l-1,ti) + 1;
            
            end
        end % loop over streams
        

    end % loop over regions
    
    end % loop over times

    % divide by the number of nodes used to calculate
    % this averages nodes that lie on the boundary of multiple regions
    uDm = uDm./divby; uDx = uDx./divby; uDy = uDy./divby; 
    if ~results.IsTwoD; uDz = uDz./divby; end
    uCm = uCm./divby; uCx = uCx./divby; uCy = uCy./divby; 
    if ~results.IsTwoD; uCz = uCz./divby; end
    Rem = Rem./divby; Rex = Rex./divby; Rey = Rey./divby; 
    if ~results.IsTwoD; Rez = Rez./divby; end
    sH = sH./divby;
    sR = sR./divby;

    % check for NaN results       
    uDm(isnan(uDm)) = 0; uDx(isnan(uDx)) = 0; uDy(isnan(uDy)) = 0;
    uCm(isnan(uCm)) = 0; uCx(isnan(uCx)) = 0; uCy(isnan(uCy)) = 0;
    Rem(isnan(Rem)) = 0; Rex(isnan(Rex)) = 0; Rey(isnan(Rey)) = 0;
    sH(isnan(sH)) = 0;
    sR(isnan(sR)) = 0;
    sR(isinf(sR)) = 0; % need this as well
    if ~results.IsTwoD
        uDz(isnan(uDz)) = 0;
        uCz(isnan(uCz)) = 0;
        Rez(isnan(Rez)) = 0;
    end

    for l = 1:l1
        
        results.DarcyFlux(:,l,:) = uDm(:,l,:);
        results.XDarcyFlux(:,l,:) = uDx(:,l,:);
        results.YDarcyFlux(:,l,:) = uDy(:,l,:);
        results.ChannelVelocity(:,l,:) = uCm(:,l,:);
        results.XChannelVelocity(:,l,:) = uCx(:,l,:);
        results.YChannelVelocity(:,l,:) = uCy(:,l,:);
        results.Reynolds(:,l,:) = Rem(:,l,:);
        results.XReynolds(:,l,:) = Rex(:,l,:);
        results.YReynolds(:,l,:) = Rey(:,l,:);
        results.StreamVolHeating(:,l,:) = sH(:,l,:);
        results.StreamVolResistance(:,l,:) = sR(:,l,:);
        
        if ~results.IsTwoD
            results.ZDarcyFlux(:,l,:) = uDz(:,l,:);
            results.ZChannelVelocity(:,l,:) = uCz(:,l,:);
            results.ZReynolds(:,l,:) = Rez(:,l,:);
        end

    end
        
    
