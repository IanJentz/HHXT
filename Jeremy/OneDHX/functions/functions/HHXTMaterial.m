function [varargout] = HHXTMaterial(model,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n1 = model.PDESystemSize;
streams = (n1-1)/2;
matlDefs = model.MaterialProperties.MaterialAssignments;
[geomDim,~] = size(model.Mesh.Nodes);

if isempty(varargin)
    varargout{1} = model.MateriallMap;
    return
end



arg = 1;
while arg <= length(varargin)
    switch varargin{arg}
        case {'BuildMaterialMap'}
            % create a map between region/stream and material propertie location
            [~,nMatl] = size(model.MaterialProperties.MaterialAssignments);
            rows = (streams+1); 
            switch geomDim
                case 3
                    cols = model.Geometry.NumCells;
                case 2
                    cols = model.Geometry.NumFaces;
            end
            matlMap = zeros(rows,cols);
            for j = 1:cols
                for i = 1:rows
                    matlInd = 0;
                    for k = 1:nMatl
                        if (model.MaterialProperties.MaterialAssignments(1,k).StreamID == (i-1)) && (model.MaterialProperties.MaterialAssignments(1,k).RegionID == j)
                            matlInd = k;
                        end
                    end
                    matlMap(i,j) = matlInd;
                end
            end
            if any(any(matlMap==0))
%                 warning('not all Regions have materials defined for all streams')
            end
            model.MaterialMap = matlMap;
            varargout{1} = matlMap;
            arg = arg+1;
            
        case {'SetAtConditions'}
            if isempty(varargin{arg+1})
                error('not enough inputs')
            end
            matlMap = model.MaterialMap;
            % based on inputs determine u0
            setDofs = false;
            setRe = false;
            u0 = zeros(n1,1); pos = []; time = [];
            switch varargin{arg+1}
                case {'pos,time,temp,press','temp,press'}
                    setDofs = true;
                    switch varargin{arg+1}
                        case 'pos,time,temp,press'
                            try
                            pos = varargin{arg+2};
                            time = varargin{arg+3};
                            temp = varargin{arg+4};
                            press = varargin{arg+5};
                            catch
                                error('not enough inputs')
                            end
                            arg = arg+6;
                        case 'temp,press'
                            pos = [0,0,0];
                            time = 0;
                            try
                            temp = varargin{arg+2};
                            press = varargin{arg+3};
                            catch
                                error('not enough inputs')
                            end
                            arg = arg+4;
                    end
                    if length(temp)==(streams+1)
                        u0(1) = temp(1);
                        for i = 1:streams
                           u0(2*i+1) = temp(i+1); 
                        end
                    elseif length(temp)==1
                        u0(1) = temp;
                        for i = 1:streams
                           u0(2*i+1) = temp; 
                        end
                    else
                        error('incorrect size for temp')
                    end
                    if length(press)==streams
                        for i = 1:streams
                           u0(2*i) = press(i); 
                        end
                    elseif length(temp)==1
                        for i = 1:streams
                           u0(2*i) = press; 
                        end
                    else
                        error('incorrect size for temp')
                    end
                    
                case {'pos,time,u0','u0'}
                    setDofs = true;
                    switch varargin{arg+1}
                        case 'pos,time,u0'
                            try
                            pos = varargin{arg+2};
                            time = varargin{arg+3};
                            u0 = varargin{arg+4};
                            catch
                                error('not enough inputs')
                            end
                            arg = arg+5;
                        case 'u0'
                            try
                            pos = [0;0;0];
                            time = 0;
                            u0 = varargin{arg+2};
                            catch
                                error('not enough inputs')
                            end
                            arg = arg+3;
                    end
                    
                case {'pos,time,Re','Re'}
                    setRe = true;
                    switch varargin{arg+1}
                        case 'pos,time,Re'
                            try
                            pos = varargin{arg+2};
                            time = varargin{arg+3};
                            Re = varargin{arg+4};
                            catch
                                error('not enough inputs')
                            end
                            arg = arg+5;
                        case 'Re'
                            pos = [0;0;0];
                            time = 0;
                            try
                            Re = varargin{arg+2};
                            catch
                                error('not enough inputs')
                            end
                            arg = arg +3;
                    end
                    
                    if length(Re)==streams
                        for i = 1:streams
                           u0(2*i) = Re(i); 
                        end
                    elseif length(Re)==1
                        for i = 1:streams
                           u0(2*i) = Re; 
                        end
                    else
                        error('incorrect size for temp')
                    end
                      
            end % switch varargin{arg+1}
            
            matlDefs = matlDefs;
            
            if setDofs==true
            %  set materials to constants at u0                    
                [rows,cols] = size(matlMap);
                for j = 1:cols
                    % for the core solid
                    temp = u0(1);
                    phi = matlDefs(matlMap(1,j)).PhaseFraction;
                    if isa(phi,'function_handle')
                        matlDefs(matlMap(1,j)).PhaseFraction = phi(pos,time); % solid materials can be position and temperature dependent
                    end
                    k = matlDefs(matlMap(1,j)).ThermalConductivity;
                    if isa(k,'function_handle')
                        matlDefs(matlMap(1,j)).ThermalConductivity = k(pos,time,temp); % solid materials can be position and temperature dependent
                    end
                    rho = matlDefs(matlMap(1,j)).MassDensity;
                    if isa(rho,'function_handle')
                        matlDefs(matlMap(1,j)).MassDensity = rho(pos,time,temp); % solid materials can be position and temperature dependent
                    end
                    cp = matlDefs(matlMap(1,j)).SpecificHeat;
                    if isa(cp,'function_handle')
                        matlDefs(matlMap(1,j)).SpecificHeat = cp(pos,time,temp); % solid materials can be position and temperature dependent
                    end
                    % for the streams
                    for i = 2:rows
                        
                    if matlMap(i,j) ~= 0
                        
                        temp = u0( 2*(i-1)+1 );
                        press = u0( 2*(i-1) );
                        phi = matlDefs(matlMap(i,j)).PhaseFraction;
                        if isa(phi,'function_handle')
                            matlDefs(matlMap(i,j)).PhaseFraction = phi(pos,time);
                        end
                        k = matlDefs(matlMap(i,j)).ThermalConductivity;
                        if isa(k,'function_handle')
                            matlDefs(matlMap(i,j)).ThermalConductivity = k(pos,time,temp,press);
                        end
                        rho = matlDefs(matlMap(i,j)).MassDensity;
                        if isa(rho,'function_handle')
                            matlDefs(matlMap(i,j)).MassDensity = rho(pos,time,temp,press);
                        end
                        cp = matlDefs(matlMap(i,j)).SpecificHeat;
                        if isa(cp,'function_handle')
                            matlDefs(matlMap(i,j)).SpecificHeat = cp(pos,time,temp,press);
                        end
                        h = matlDefs(matlMap(i,j)).Enthalpy;
                        if isa(h,'function_handle')
                            matlDefs(matlMap(i,j)).Enthalpy = h(pos,time,temp,press);
                        end
                        mu = matlDefs(matlMap(i,j)).Viscosity;
                        if isa(mu,'function_handle')
                            matlDefs(matlMap(i,j)).Viscosity = mu(pos,time,temp,press);
                        end
                        Dh = matlDefs(matlMap(i,j)).HydraulicDiameter;
                        if isa(Dh,'function_handle')
                            matlDefs(matlMap(i,j)).HydraulicDiameter = Dh(pos,time,temp,press);
                        end
                        
                    end
                    end % loop over rows (streams)
                end % loop over cols
            end
            
            if setRe==true
                [rows,cols] = size(matlMap);
                for j = 1:cols
                    % for the streams
                    for i = 2:rows
                        
                    if matlMap(i,j) ~= 0
                        
                        Re = u0( 2*(i-1) );
                        f = matlDefs(matlMap(i,j)).fDarcy;
                        if isa(f,'function_handle')
                            matlDefs(matlMap(i,j)).fDarcy = f(Re*ones(3,1));
                        end
                        jCol = matlDefs(matlMap(i,j)).jColburn;
                        if isa(jCol,'function_handle')
                            matlDefs(matlMap(i,j)).jColburn= jCol(pos,time,Re);
                        end
                        R = matlDefs(matlMap(i,j)).Resistance;
                        if isa(R,'function_handle')
                            matlDefs(matlMap(i,j)).Resistance = R(pos,time,Re);
                        end
                        kp = matlDefs(matlMap(i,j)).Permeability;
                        if isa(kp,'function_handle')
                            matlDefs(matlMap(i,j)).Permeability = kp(pos,time,Re);
                        end
                        
                    end
                    end % loop over rows (streams)
                end % loop over cols
                            
            end
            
            varargout{1} = matlDefs;
            model.MaterialProperties.MaterialAssignments = matlDefs;
            
        case {'SetMaterialDefinitions'}
            if ~isa(varargin{arg+1},'HHXTMaterialAssignment')
                error('input is not of class HHXTMaterialAssignment')
            end
            model.MaterialProperties.MaterialAssignments = varargin{arg+1};
            varargin{end+1} = 'BuildMaterialMap';
            arg = arg+2;
            
        otherwise
            error('unrecognized input')
            
            
    end % switch varargin{arg}
end % loop over varargin



end

