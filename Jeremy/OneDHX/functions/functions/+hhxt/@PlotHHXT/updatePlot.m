function hp = updatePlot(obj,plottype,varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

figure(obj.fighandle);
if nargin >= 3
obj = obj.setOptions(varargin);
end

n1 =  obj.model.PDESystemSize;
[geomDims,nNodes] = size(obj.model.Mesh.Nodes);



switch obj.model.AnalysisType
    
    case 'steadystate'
        
        switch geomDims
            %% Steady State 3D
            case 3
                switch plottype
                    
                    
                end % switch plottype
                
            %% Steady State 2D
            case 2
                switch plottype
                    
                    case 'NodalSolutions'
                        u = obj.results.NodalSolution;
%                         u = regionset(u,obj.options.nodes);
                        u = hideUndefined(u,obj.j_nostr);
                        for i = 1:n1
                            subplot(n1,1,i);
                            if ~isempty(obj.options.addArgs)
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,i),obj.options.addArgs{1,:});
                            contourLines(hp,obj.options.addArgs)
                            else
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,i));    
                            end
                            if i == 1
                                title('Body : $T$ $[C]$','interpreter','latex');
                                setCLim(obj.options.templimits);
                            elseif ~mod(i,2)
                                title(['Fluid ',num2str(i/2),' : $P$ $[Pa]$'],'interpreter','latex');
                                setCLim(obj.options.presslimits);
                            else
                                title(['Fluid ',num2str((i-1)/2),' : $T$ $[C]$'],'interpreter','latex');
                                setCLim(obj.options.templimits);
                            end
                            overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                            setView([]);
                        end
                        
                    case 'NodalGradients'
                        dudx = obj.results.XGradients;
                        dudy = obj.results.YGradients;
                        dudx = hideUndefined(dudx,obj.j_nostr);
                        dudy = hideUndefined(dudy,obj.j_nostr);
                        dud = [dudx;dudy];
                        for i = 1:n1
                            cclimits = findCLim(obj.options.climits,dud(:,i));
                            subplot(n1,2,2*i-1);
                            if ~isempty(obj.options.addArgs)
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',dudx(:,i),obj.options.addArgs{1,:});
                            contourLines(hp,obj.options.addArgs)
                            else
                            pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',dudx(:,i));  
                            end
                            if i == 1
                                title('Body : $\frac{dT}{dx}$ $[\frac{K}{m}]$','interpreter','latex');
                            elseif ~mod(i,2)
                                title(['Fluid ',num2str(i/2),' : $\frac{dP}{dx}$ $[\frac{Pa}{m}]$'],'interpreter','latex');
                            else
                                title(['Fluid ',num2str((i-1)/2),' : $\frac{dT}{dx}$ $[\frac{K}{m}]$'],'interpreter','latex');
                            end
                            setCLim(cclimits);
                            overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                            setView([]);
                            subplot(n1,2,2*i);
                            if ~isempty(obj.options.addArgs)
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',dudy(:,i),obj.options.addArgs{1,:});
                            contourLines(hp,obj.options.addArgs)
                            else
                            pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',dudy(:,i));  
                            end
                            if i == 1
                                title('Body : $\frac{dT}{dy}$ $[\frac{K}{m}]$','interpreter','latex');
                            elseif ~mod(i,2)
                                title(['Fluid ',num2str(i/2),' : $\frac{dP}{dy}$ $[\frac{Pa}{m}]$'],'interpreter','latex');
                            else
                                title(['Fluid ',num2str((i-1)/2),' : $\frac{dT}{dy}$ $[\frac{K}{m}]$'],'interpreter','latex');
                            end
                            setCLim(cclimits);
                            overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                            setView([]);
                        end
                        
                        case 'NodalGradientsMagnitude'
                            dudx = obj.results.XGradients;
                            dudy = obj.results.YGradients;
                            dudx = hideUndefined(dudx,obj.j_nostr);
                            dudy = hideUndefined(dudy,obj.j_nostr);
                            dud = sqrt(dudx.^2+dudy.^2);
                            for i = 1:n1
                                cclimits = findCLim(obj.options.climits,dud(i,:));
                                subplot(n1,1,i);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',dud(:,i));
                                if i == 1
                                    title('Body : $\left |\nabla T  \right |$ $[\frac{K}{m}]$','interpreter','latex');
                                elseif ~mod(i,2)
                                    title(['Fluid ',num2str(i/2),' : $\left |\nabla P  \right |$ $[\frac{Pa}{m}]$'],'interpreter','latex');
                                else
                                    title(['Fluid ',num2str((i-1)/2),' : $\left |\nabla T  \right |$ $[\frac{K}{m}]$'],'interpreter','latex');
                                end
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                        case 'DarcyVelocityMagnitude'
                            delete(obj.fighandle.Children)
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                uu = obj.results.DarcyFlux(:,i);
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),1,i);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                title(['Fluid ',num2str(i),' : $\left | u_{d} \right |$ $[\frac{m}{s}]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                        case 'ChannelVelocity'
                            delete(obj.fighandle.Children)
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                ux = obj.results.XChannelVelocity(:,i);
                                uy = obj.results.YChannelVelocity(:,i);
                                ux = hideUndefined(ux,obj.j_nostr,ind(i));
                                uy = hideUndefined(uy,obj.j_nostr,ind(i));
                                uu = [ux;uy];
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),2,2*(i-1)+1);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',ux);
                                title(['Fluid ',num2str(i),' : $u_{f,x}$ $[\frac{m}{s}]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                                subplot(length(ind),2,2*(i-1)+2);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uy);
                                title(['Fluid ',num2str(i),' : $u_{f,y}$ $[\frac{m}{s}]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                        case 'ChannelVelocityMagnitude'
                            delete(obj.fighandle.Children)
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                uu = obj.results.ChannelVelocity(:,i);
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),1,i);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                title(['Fluid ',num2str(i),' : $\left | u_{f} \right |$ $[\frac{m}{s}]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end


                        case 'ChannelVelocityVector'
                            delete(obj.fighandle.Children)
                            node_pos = obj.model.Mesh.Nodes';
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                ux = obj.results.XChannelVelocity(:,i);
                                uy = obj.results.YChannelVelocity(:,i);
                                ux = hideUndefined(ux,obj.j_nostr,ind(i));
                                uy = hideUndefined(uy,obj.j_nostr,ind(i));
                                uu = obj.results.ChannelVelocity(:,i);
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                subplot(length(ind),1,i);
                                if ~isempty(obj.options.addArgs)
                                quiver(node_pos(:,1),node_pos(:,2),ux,uy,obj.options.addArgs{:})
                                else
                                quiver(node_pos(:,1),node_pos(:,2),ux,uy)    
                                end
%                                     pdeplot(obj.model,'FlowData',[ux,uy]);
                                title(['Fluid ',num2str(i),' : $\left | u_{f} \right |$ $[\frac{m}{s}]$'],'interpreter','latex');
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                            case 'Reynolds+ChannelVelocityVector'
                            delete(obj.fighandle.Children)
                            node_pos = obj.model.Mesh.Nodes';
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                ux = obj.results.XChannelVelocity(:,i);
                                uy = obj.results.YChannelVelocity(:,i);
                                ux = hideUndefined(ux,obj.j_nostr,ind(i));
                                uy = hideUndefined(uy,obj.j_nostr,ind(i));
                                uu = obj.results.Reynolds(:,i);
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),1,i);
                                
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                setCLim(cclimits);
                                hold on
                                if ~isempty(obj.options.addArgs)
                                quiver(node_pos(:,1),node_pos(:,2),ux,uy,obj.options.addArgs{:})
                                else
                                quiver(node_pos(:,1),node_pos(:,2),ux,uy)    
                                end
%                                     pdeplot(obj.model,'FlowData',[ux,uy]);
                                hold off
                                title(['Fluid ',num2str(i),' : $\left | Re \right |$ $[-]$'],'interpreter','latex');
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                            case 'Press+ChannelVelocityVector'
                            delete(obj.fighandle.Children)
                            node_pos = obj.model.Mesh.Nodes';
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                ux = obj.results.XChannelVelocity(:,i);
                                uy = obj.results.YChannelVelocity(:,i);
%                                 logminx = min(log(abs(ux(ux~=0))));
%                                 logminy = min(log(abs(uy(uy~=0)))) ;
%                                 ux = (log(abs(ux))-logminx).*sign(ux);
%                                 uy = (log(abs(uy))-logminy).*sign(uy);
                                ux = sqrt(abs(ux)).*sign(ux);
                                uy = sqrt(abs(uy)).*sign(uy);
                                ux = hideUndefined(ux,obj.j_nostr,ind(i));
                                uy = hideUndefined(uy,obj.j_nostr,ind(i));
                                uu = obj.results.NodalSolution(:,(i*2));
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),1,i);
                                
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                setCLim(cclimits);
                                hold on
                                if ~isempty(obj.options.addArgs)
                                quiver(node_pos(:,1),node_pos(:,2),ux,uy,obj.options.addArgs{:})
                                else
                                quiver(node_pos(:,1),node_pos(:,2),ux,uy)    
                                end
%                                     pdeplot(obj.model,'FlowData',[ux,uy]);
                                hold off
                                title(['Fluid ',num2str(i),' :',newline,'vectors of $\mathbf{v}$ \& color map of $P$ $[Pa]$'],'interpreter','latex');
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                                
                            case 'Reynolds'
                            delete(obj.fighandle.Children)
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                ux = obj.results.XReynolds(:,i);
                                uy = obj.results.YReynolds(:,i);
                                ux = hideUndefined(ux,obj.j_nostr,ind(i));
                                uy = hideUndefined(uy,obj.j_nostr,ind(i));
                                uu = [ux;uy];
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),2,2*(i-1)+1);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',ux);
                                title(['Fluid ',num2str(i),': $Re_{x}$ $[-]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                                subplot(length(ind),2,2*(i-1)+2);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uy);
                                title(['Fluid ',num2str(i),' : $Re_{y}$ $[-]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                            case 'ReynoldsMagnitude'
                                delete(obj.fighandle.Children)
                                ind = 1:n1;
                                ind = ind(~mod(ind,2));
                                for i = 1:length(ind)
                                    uu = obj.results.Reynolds(:,i);
                                    uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                    cclimits = findCLim(obj.options.climits,uu);
                                    subplot(length(ind),1,i);
                                    pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                    title(['Fluid ',num2str(i),' : $\left | Re \right |$ $[-]$'],'interpreter','latex');
                                    setCLim(cclimits);
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                            
                            
                            case 'ReynoldsVector'
                                delete(obj.fighandle.Children)
                                ind = 1:n1;
                                ind = ind(~mod(ind,2));
                                for i = 1:length(ind)
                                    ux = obj.results.XReynolds(:,i);
                                    uy = obj.results.YReynolds(:,i);
                                    ux = hideUndefined(ux,obj.j_nostr,ind(i));
                                    uy = hideUndefined(uy,obj.j_nostr,ind(i));
                                    uu = obj.results.ChannelVelocity(:,i);
                                    uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                    subplot(length(ind),1,i);
                                    pdeplot(obj.model,'FlowData',[ux,uy]);
                                    title(['Fluid ',num2str(i),' : $Re$ $[-]$'],'interpreter','latex');
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                                
                            case 'Temperatures'
                                delete(obj.fighandle.Children)
                                u = obj.results.NodalSolution;
        %                         u = regionset(u,obj.options.nodes);
                                u = hideUndefined(u,obj.j_nostr);
                                tind = 1:2:n1;
                                for i = 1:length(tind)
                                    subplot((n1+1)/2,1,i);
                                    if ~isempty(obj.options.addArgs)
                                    hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)),obj.options.addArgs{1,:});
                                    contourLines(hp,obj.options.addArgs)
                                    else
                                    hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)));    
                                    end
                                    if i == 1
                                        title('Body : $T$ $[C]$','interpreter','latex');
                                        setCLim(obj.options.templimits);
                                    else
                                        title(['Fluid ',num2str((i-1)),' : $T$ $[C]$'],'interpreter','latex');
                                        setCLim(obj.options.templimits);
                                    end
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                                
                            case 'FluidTemperatures'
                                delete(obj.fighandle.Children)
                                u = obj.results.NodalSolution;
        %                         u = regionset(u,obj.options.nodes);
                                u = hideUndefined(u,obj.j_nostr);
                                tind = 1:2:n1;
                                for i = 2:length(tind)
                                    subplot(((n1+1)/2)-1,1,i-1);
                                    if ~isempty(obj.options.addArgs)
                                    hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)),obj.options.addArgs{1,:});
                                    contourLines(hp,obj.options.addArgs)
                                    else
                                    hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)));    
                                    end
                                    if i == 1
                                        title('Body : $T$ $[C]$','interpreter','latex');
                                        setCLim(obj.options.templimits);
                                    else
                                        title(['Fluid ',num2str((i-1)),' : $T$ $[C]$'],'interpreter','latex');
                                        setCLim(obj.options.templimits);
                                    end
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                                
                            case 'BodyTemperature'
                                delete(obj.fighandle.Children)
                                u = obj.results.NodalSolution;
        %                         u = regionset(u,obj.options.nodes);
                                u = hideUndefined(u,obj.j_nostr);
                                tind = 1:2:n1;
                                for i = 1
                                    if ~isempty(obj.options.addArgs)
                                    hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)),obj.options.addArgs{1,:});
                                    contourLines(hp,obj.options.addArgs)
                                    else
                                    hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)));    
                                    end

                                        title('Body : $T$ $[C]$','interpreter','latex');
                                        setCLim(obj.options.templimits);
            
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                                
                            case 'FluidHeating'
                                delete(obj.fighandle.Children)
                                ind = 1:n1;
                                ind = ind(~mod(ind,2));
                                for i = 1:length(ind)
                                    uu = obj.results.StreamVolHeating(:,i);
                                    uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                    sum(uu(~isnan(uu)));
                                    cclimits = findCLim(obj.options.climits,uu);
                                    subplot(length(ind),1,i);
                                    pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                    title(['Fluid ',num2str(i),' : $\dot{Q}$ $[\frac{W}{m^3}]$'],'interpreter','latex');
                                    setCLim(cclimits);
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                                
                            case 'Resistance'
                                delete(obj.fighandle.Children)
                                ind = 1:n1;
                                ind = ind(~mod(ind,2));
                                for i = 1:length(ind)
                                    uu = obj.results.StreamVolResistance(:,i);
                                    uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                    sum(uu(~isnan(uu)))
                                    cclimits = findCLim(obj.options.climits,uu);
                                    subplot(length(ind),1,i);
                                    pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                    title(['Fluid ',num2str(i),' : $R_{V}$ $[\frac{Km^3}{W}]$'],'interpreter','latex');
                                    setCLim(cclimits);
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                    
                    
                end % switch plottype
                
        end % switch geomDIms
        
        
        
    case 'transient'
        
        t_ind = obj.options.time;
        if length(t_ind) > 1
            t_vind = t_ind;
            t_ind = t_ind(1);
        end
        
        switch geomDims
            %% Transient 3D
            case 3
                switch plottype
                    
                    
                end % switch plottype
                
            %% Transient 2D
            case 2
                switch plottype
                    
                    case 'NodalSolutions'
                        u = obj.results.NodalSolution(:,:,t_ind);
%                         u = regionset(u,obj.options.nodes);
                        u = hideUndefined(u,obj.j_nostr);
                        for i = 1:n1
                            subplot(n1,1,i);
                            if ~isempty(obj.options.addArgs)
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,i),obj.options.addArgs{1,:});
                            contourLines(hp,obj.options.addArgs)
                            else
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,i));    
                            end
                            if i == 1
                                title('Core: $T$ $[C]$','interpreter','latex');
                                setCLim(obj.options.templimits);
                            elseif ~mod(i,2)
                                title(['Fluid ',num2str(i/2),': $P$ $[Pa]$'],'interpreter','latex');
                                setCLim(obj.options.presslimits);
                            else
                                title(['Fluid ',num2str((i-1)/2),': $T$ $[C]$'],'interpreter','latex');
                                setCLim(obj.options.templimits);
                            end
                            overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                            setView([]);
                        end
                        
                    case 'Temperatures'
                        delete(obj.fighandle.Children)
                        u = obj.results.NodalSolution(:,:,t_ind);
%                         u = regionset(u,obj.options.nodes);
                        u = hideUndefined(u,obj.j_nostr);
                        tind = 1:2:n1;
                        for i = 1:length(tind)
                            subplot((n1+1)/2,1,i);
                            if ~isempty(obj.options.addArgs)
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)),obj.options.addArgs{1,:});
                            contourLines(hp,obj.options.addArgs)
                            else
                            hp = pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',u(:,tind(i)));    
                            end
                            if i == 1
                                title('Core: $T$ $[C]$','interpreter','latex');
                                setCLim(obj.options.templimits);
                            else
                                title(['Fluid ',num2str((i-1)),': $T$ $[C]$'],'interpreter','latex');
                                setCLim(obj.options.templimits);
                            end
                            overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                            setView([]);
                        end
                        
                    case 'ReynoldsMagnitude'
                                delete(obj.fighandle.Children)
                                ind = 1:n1;
                                ind = ind(~mod(ind,2));
                                for i = 1:length(ind)
                                    uu = obj.results.Reynolds(:,i,t_ind);
                                    uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                    cclimits = findCLim(obj.options.climits,uu);
                                    subplot(length(ind),1,i);
                                    pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                    title(['Fluid ',num2str(i),': $\left | Re \right |$ $[-]$'],'interpreter','latex');
                                    setCLim(cclimits);
                                    overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                    setView([]);
                                end
                                
                    case 'DarcyVelocityMagnitude'
                            delete(obj.fighandle.Children)
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                uu = obj.results.DarcyFlux(:,i,t_ind);
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),1,i);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                title(['Fluid ',num2str(i),': $\left | u_{d} \right |$ $[\frac{m}{s}]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                            
                    case 'FluidHeating'
                            delete(obj.fighandle.Children)
                            ind = 1:n1;
                            ind = ind(~mod(ind,2));
                            for i = 1:length(ind)
                                uu = obj.results.StreamVolHeating(:,i,t_ind);
                                uu = hideUndefined(uu,obj.j_nostr,ind(i));
                                cclimits = findCLim(obj.options.climits,uu);
                                subplot(length(ind),1,i);
                                pdeplot(obj.model,'ColorMap',obj.options.clrmap,'XYData',uu);
                                title(['Fluid ',num2str(i),': $\dot{Q}$ $[\frac{W}{m}]$'],'interpreter','latex');
                                setCLim(cclimits);
                                overlayPlots2D(obj.model,obj.options.wireframe,obj.options.meshplot);
                                setView([]);
                            end
                    
                end % switch plottype
                
        end % switch geomDIms
        
end % switch AnalysisType





%% Internal Functions
    function v = setView(vin)
        if isempty(vin)
        view([0,0]);
%         view(90,-90);
        view(0,90);
        else
            view(vin);
        end
        [ax,el] = view;
        v = [ax,el];
    end

    function setCLim(limits)
        if ~isempty(limits)
            caxis(limits);
        end
    end

    function climits = findCLim(limits,u)
        if isempty(limits)
            climits = [];
            return
        end
        if ~ischar(limits)
            climits = limits;
            return
        end
        switch limits
            case 'magnitude'
                climits = [min(u),max(u)];
            case 'centralmagnitude'
                u = max(abs(min(u)),max(u));
                climits = [-u,u];
            otherwise
                climits = limits;
        end
    end

    function setColorMap(colormapopts)
        if isempty(colormapopts)
            
        else
            switch colormapopts
                case 'single'
                    colorbar('off');
                case 'west'
                    colorbar('west')
                case 'on'
                    colorbar()
            end
        end
    end

    function overlayPlots2D(model,wireframe,meshplot)
        
        
        hold on
        if (wireframe == true) || (meshplot == true)
            h = pdegplot(model);
            h.Color = [0,0,0];
        end
        if meshplot == true
            hold on
            h = pdeplot(model);
            line = h.findobj;
            line(1).Color = [0,0,0,0.25];
            line(2).Color = [0,0,0];
        end
        
        hold off
        
    end

    function contourLines(h,addArgs)
       cont = false;
       for iarg = 1:length(addArgs)
           switch addArgs{iarg}
               case 'Levels'
                   nlvls = addArgs{iarg+1};
               case 'Contour'
                   cont = true;
           end
       end
       
       if cont == true
           line = h.findobj;
           for l = 2:nlvls+1
               line(l).Color = [1,1,1];
           end
       end
        
    end

    function uu = regionset(u,nodes)
       
        if ~isempty(nodes)
            [nr,nc] = size(u);
            uu = ones(nr,nc)*NaN;
            uu(nodes,:) = u(nodes,:);
        else
            uu = u;
        end
        
    end

    function u = hideUndefined(u,j_nostr,col)
        [~,cols] = size(u);
        if nargin < 3
            for col = 1:cols
            if col > 1; u(j_nostr(:,ceil((col-1)/2)),col) = NaN; end
            end
        else
            u(j_nostr(:,ceil((col-1)/2)),1) = NaN;
        end
    end
end

