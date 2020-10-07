classdef StreamBoundaryCondition < handle & matlab.mixin.internal.CompactDisplay
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        StreamID; %stream index - 0 for solid, 1, 2, etc. for fluid streams
        ThermalTypes; %type of thermal boundary - 'Temperature', 'HeatFlux', 'HeatFlow'
        HydraulicTypes; %type of hydraulic boundary - 'Pressure', 'MassFlux', 'MassFlow'
        ThermalVals; %value of the thermal BC
        HydraulicVals; %value of the hydraulic BC
        FlowDir; %direction of flow - 'inlet', or 'outlet'
        RegionType; %'Face' or 'Cell'
        RegionID; % region ID of the boundary
        BCarea; % area in m2 of the boundary
        BCType; % matlab BC type - 'dirichlet', 'neumann', 'mixed'
        nodes; % nodes of the boundary
        elems; % elements of the boundary
       
    end
    
    methods(Hidden=true)
        function self = StreamBoundaryCondition(bcr,varargin)
            self.RecordOwner = bcr;
            pdem = bcr.ParentPdemodel;
            self.ParentPdemodel = pdem;
            parser = inputParser;
            parser.addParameter('SystemSize', []);
            parser.addParameter('Edge', []);
            parser.addParameter('Face', []);
            parser.addParameter('Stream', []);
            parser.addParameter('Direction', []);
            parser.addParameter('Temperature', []);
            parser.addParameter('HeatFlux', []);
            parser.addParameter('HeatFlow', []);
            parser.addParameter('Pressure', []);
            parser.addParameter('MassFlux', []);
            parser.addParameter('MassFlow', []);
            parser.parse(varargin{:});
        
            if ~isempty(parser.Results.Face)
                self.RegionType = 'Face';
                self.RegionID = parser.Results.Face;
                self.nodes = pdem.Mesh.findNodes('region','Face',self.RegionID);
                self.elems = pdem.Mesh.findElements('region','Face',self.RegionID);
                self.BCarea = area(pdem.Mesh,self.elems);
            elseif ~isempty(parser.Results.Edge)
                self.RegionType = 'Edge';
                self.RegionID = parser.Results.Edge;
                self.nodes = pdem.Mesh.findNodes('region','Edge',self.RegionID);

                [pm,em,~] = pdem.Mesh.meshToPet;
                self.elems = find(em(5,:)==self.RegionID);
                sem = em(:,self.elems);
                elen = 0;
                for di = 1:size(sem,2)
                    elen = elen + norm(pm(:,sem(2,di))-pm(:,sem(1,di)));
                end
                self.BCarea = pdem.Thickness*elen;
                
            end
            
            dirmult = 1;
            if ~isempty(parser.Results.Stream)
                str = parser.Results.Stream; 
                if ~isempty(parser.Results.Direction)
                    self.FlowDir = parser.Results.Direction;
                    if strcmpi(parser.Results.Direction,'outlet')
                        dirmult = -1;
                        self.FlowDir = 'outlet';
                    elseif strcmpi(parser.Results.Direction,'inlet')
                        self.FlowDir = 'inlet';
                    else
                        self.FlowDir = [];
                    end
                end
            else
                str = 0;
            end
            self.StreamID = str;
            
            
            ni = ~isempty(parser.Results.Temperature)...
                + ~isempty(parser.Results.HeatFlux) + ~isempty(parser.Results.HeatFlow);
            nj = ~isempty(parser.Results.Pressure) + ~isempty(parser.Results.MassFlux)...
                + ~isempty(parser.Results.MassFlow);
            if ni ~= 0
            self.ThermalTypes = cell(1,ni);
            self.ThermalVals = cell(1,ni);
            end
            if nj ~= 0
            self.HydraulicTypes = cell(1,nj);
            self.HydraulicVals = cell(1,nj);
            end
            
            % parse the stream BCs           
            di = 1; ni = 1;
            if ~isempty(parser.Results.Temperature)
                self.ThermalTypes{di} = 'Temperature';
                self.ThermalVals{di} = parser.Results.Temperature;
                di = di+1;
            end
            if ~isempty(parser.Results.HeatFlux)
                self.ThermalTypes{di} = 'HeatFlux';
                self.ThermalVals{di} = parser.Results.HeatFlux;
                di = di+1;
            end
            if ~isempty(parser.Results.HeatFlow)
                self.ThermalTypes{di} = 'HeatFlow';
                self.ThermalVals{di} = parser.Results.HeatFlow;
                di = di+1;
            end
            if ~isempty(parser.Results.Pressure)
                self.HydraulicTypes{ni} = 'Pressure';
                self.HydraulicVals{ni} = parser.Results.Pressure;
                ni = ni+1;
            end
            if ~isempty(parser.Results.MassFlux)
                self.HydraulicTypes{ni} = 'MassFlux';
                self.HydraulicVals{ni} = parser.Results.MassFlux;
                ni = ni+1;
            end
            if ~isempty(parser.Results.MassFlow)
                self.HydraulicTypes{ni} = 'MassFlow';
                self.HydraulicVals{ni} = parser.Results.MassFlow;
                ni = ni+1;       
            end
                      
            
        
        end
        
       
    end
    
    methods
        function T_mean = meanTemperature(self)
            %MEANTEMPERATURE Element mean temperature at the boundary
            n = self.StreamID*2+1;
            T_mean = sum(self.ParentPdemodel.Runtime.ugm(n,self.nodes));
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    
    properties (Hidden = true, Access=?hhxt.HHXTModel)
%     (Hidden = true, SetAccess='private')
        ParentPdemodel;
    end  
    
    properties (Hidden = true, Access='private')
        RecordOwner;
    end  
  
end

