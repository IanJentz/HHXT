classdef PlotHHXT
    %PLOTHHXT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fighandle
        model
        results
        j_nostr
        options
    end
    
    methods
        function obj = PlotHHXT(filename,fighandle)
            %PLOTHHXT Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 1
                obj.fighandle = figure('Position',[100 100 1400 800]);
            else
                obj.fighandle = fighandle;
            end
            load(filename,'model','results');
            obj.model = model;
            
            n1 =  model.PDESystemSize;
            [~,nNodes] = size(model.Mesh.Nodes);
            n_str = (model.PDESystemSize-1)/2;
            
            obj.results = results;
            obj = obj.resetOptions();
            
            % make a boolean array indicating nodes that are not part of streams
            obj.j_nostr = true(nNodes,n_str);
            for str = 1:n_str
                m = model.MaterialMap(str+1,(model.MaterialMap(str+1,:)~=0));
                for i = 1:length(m)
                    matl = model.MaterialProperties.MaterialAssignments(m(i));
                    j = findNodes(model.Mesh,'region',matl.RegionType,matl.RegionID);
                    obj.j_nostr(j,str) = false;
                end
            end
        end
        
        function obj = resetOptions(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.options.nodes = [];
            obj.options.templimits = [];
            obj.options.presslimits = [];
            obj.options.climits = [];
            obj.options.colormapopts = [];
            obj.options.addArgs = {};
            obj.options.clrmap = 'jet';
            obj.options.time = 0;
            obj.options.overlay = false;
            obj.options.wireframe = false;
            obj.options.meshplot = false;
            obj.options.showUndefined = false;
        end
        
        function obj = setOptions(obj,varargin)
            
            if nargin == 1
                [geomDims,nNodes] = size(obj.model.Mesh.Nodes);
                switch geomDims
                    case 3
                    fprintf('''Cell'',1\n')
                    case 2
                    fprintf('''Face'',1\n')
                end
                fprintf('''TempLimits'',[T_low,T_high]\n')
                fprintf('''PressLimits'',[P_low,P_high]\n')
                fprintf('''CLimits'',{[C_low,C_high],''Magnitude'',''CenteredMagnitude''}\n')
                fprintf('''Time'',time\n')
                fprintf('''Wireframe'',{''on'',''off''}\n')
                fprintf('''Mesh'',{''on'',''off''}\n')
                fprintf('''Name'',''Value'' pairs in pdeplot')
            end
            
            % when passing a cell array of arguments in it always shows up in the
            % function as a cell array within a 1x1 cell varargin
            if (isa(varargin,'cell') && length(varargin) == 1 )
                varargin = varargin{1};
            end

            if mod(length(varargin),2)
                error(['input missmatch, argument expected after ',varargin{length(nvarargin)}])
            end
            for i = 1:2:length(varargin)
                switch varargin{i}
                    case {'Cell'}
                        obj.options.nodes = obj.options.model.Mesh.findNodes('region','Cell',varargin{i+1});
                    case {'Face'}
                        obj.options.nodes = obj.options.model.Mesh.findNodes('region','Face',varargin{i+1});
                    case {'TempLimits'}
                        obj.options.templimits = varargin{i+1};
                    case {'PressLimits'}
                        obj.options.templimits = varargin{i+1};
                    case {'CLimits'}
                        if ~isa(varargin{i+1},'char')
                            obj.options.climits = varargin{i+1};
                        else
                        switch varargin{i+1}
                            case 'Magnitude'
                                obj.options.climits = 'magnitude';
                                obj.options.colormapopts = 'single';
                            case {'CenteredMagnitude','MagnitudeCentered'}
                                obj.options.climits = 'centralmagnitude';
                                obj.options.colormapopts = 'single';
                        end
                        end
                    case {'Time'}
                        obj.options.time = varargin{i+1};
                    case {'Wireframe'}
                        switch varargin{i+1}
                            case 'on'
                                obj.options.wireframe = true;
                                obj.options.overlay = true;
                            case 'off'
                                obj.options.wireframe = false;
                                obj.options.overlay = false;

                        end
                    case {'Mesh'}
                        switch varargin{i+1}
                            case 'on'
                                obj.options.meshplot = true;
                                obj.options.overlay = true;
                            case 'off'
                                obj.options.meshplot = false;
                                obj.options.overlay = false;

                        end
                    case {'pdeplotArgs'}
                        obj.options.addArgs =  varargin{i+1};

                end
            end
            
        end
        
        fig = updatePlot(obj,plottype,varargin);
            
            
        
    end
end

