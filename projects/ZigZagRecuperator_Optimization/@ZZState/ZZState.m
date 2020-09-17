classdef ZZState
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        bdir
        wdir
        ddir
        fileData
        fileFiber
        time_start
        seconds
        ntimes
        sec_diff % time_start - fiber time_start
        
        
        m_dot_H
        m_dot_C
        
        T_H_in
        T_H_out
        T_C_in
        T_C_out
        
        T_W
        pos_W
        
        Fiber
        pos_F
        T_F
        pos_TF
        pos_Map
        FilterFiber
        
        P_H_in
        DP_H
        P_C_in
        DP_C
        
        EES_sol
        EES_sol_units
        EES_arr
        EES_err
        
        debug

    end
    
    properties (Access = private)
        x_headers
        y_headers
        x_body
        y_body
    end
    
    methods
        function obj = ZZState(fileData,fileFiber,debug)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            
            if nargin == 2
                debug = false;
            end
            obj.debug = debug;
            obj.fileData = fileData;
            obj.fileFiber = fileFiber;
            obj.FilterFiber = false;
            
            obj.bdir = 'W:\PCHE\AFandZZ';
            obj.wdir = [obj.bdir,'\postProcessing'];
            
            obj = obj.setDrawingVars;
        end
        
   
        function obj = readFiber(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            load([obj.wdir,'\ZZ_Fsubsections.mat'],'subsections');
            obj = obj.readFiberFile([obj.wdir,'\',obj.fileFiber],subsections);
            
            
            %choose "shortest" fiber
            nF = length(obj.Fiber.subsecPos);
            nFp_min = Inf;
            for i = 1:nF
                nFp = size(obj.Fiber.subsecInd{i},2);
                if nFp < nFp_min
                    nFp_min = nFp;
                    dx = obj.Fiber.subsecPos{i}(2)-obj.Fiber.subsecPos{i}(1);
                end
            end
            for i = 1:nF
                obj.Fiber.subsecInd{i} = obj.Fiber.subsecInd{i}(1:nFp_min);
                obj.Fiber.subsecPos{i}(2) = obj.Fiber.subsecPos{i}(1)+dx;
            end
            obj.pos_Map = struct('x',1,'y',1);
            obj.pos_Map.x =  linspace(obj.pos_F(1,1),obj.pos_F(1,1)-dx,nFp_min);
            obj.pos_Map.y = obj.pos_F(2,:);
           
        end
    
        obj = readData(obj)
        obj = calcJ(obj,date_start,elapsed_sec)
        
        hfig = initMap(obj)
        hfig = updateMap(obj,hfig,seconds)
        v = videoMap(obj)
        hfig = initFiberPlot(obj)
        hfig = updateFiberPlot(obj,hfig,seconds)
        v = videoFiberPlot(obj)

    end
    
    methods (Access = private)
        obj = setDrawingVars(obj)
        FiberReading = readFiberFile(obj,file,subsections)
    end
    
end

