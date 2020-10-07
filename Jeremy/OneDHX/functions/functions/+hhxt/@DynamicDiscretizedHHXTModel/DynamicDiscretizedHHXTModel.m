classdef DynamicDiscretizedHHXTModel < pde.DiscretizedPDEModel
%DiscretizedPDEModel - Finite Element form of PDE
%This undocumented function may be removed in a future release.
%
%       Copyright 2015-2017 The MathWorks, Inc.
    
    properties
        
        coefstructb;
        vb;
        Ni;  % dofs that reside at nodes within undefined regions.  These are not computed
        numIgnoredEqns; % number opf dofs that are ignored ^ see above
        Bkeep;
        outflow_ind;
        
        C; %global damping matrix
        Mass; %MassMatrix
        IsSecondOrderODE;
        tspan;
        % Flags indicating whether coefficient is dependent on time or u
        vc, va, vf, vq, vg, vh, vr;
        % Flag to indicate if 'm' or 'd' coefficients depend on time.
        vd; 
        % Flags indicating whether Dirichlet BC is function of time,
        % useful in recovering full solution after 
        ht, rt;
        IsSpatialCoefficientsNonlinear;
        dudt;
        
      
    end
    
    properties(Access=private)
        % Flags indicating whether coefficients are functions of time
        ct, at, ft, qt, gt;
        dt, du
    end
    
    methods
        function obj=DynamicDiscretizedHHXTModel(thePde,p,e,t,coefstruct,coefstructb,u0,tlist,tsecondOrder)
            obj=obj@pde.DiscretizedPDEModel(thePde,p,e,t,coefstruct,u0);
            
            obj = obj.checkQforUDependence(u0,0);
            
            obj.coefstructb = coefstructb;  % add in b coefficient
            obj.qt=false; obj.gt=false; obj.ht=false; obj.rt=false; obj.ct=false;
            obj.at=false; obj.ft=false;
            obj.dt=false; obj.du=false;
            obj.IsSpatialCoefficientsNonlinear = false;
            obj.IsSecondOrderODE = tsecondOrder;
            obj.checkTlist(tlist)
            obj = obj.checkSpatialCoefsForUorTDependence(u0,tlist);
            
            obj.vc = true; % have to add this to force K matrix to be recalculated each time step
%             obj.vq = true; % have to add this to force the outflow condition to be calculated each time step
        end
    end
    
    methods(Access=protected)
        function self=checkSpatialCoefsForUorTDependence(self,u0,tlist)
            %Determine whether coefficients are functions of t
            
            t0 = min(tlist);
            if(t0 == 0)
                t0 = eps;
            end
            t1 = NaN;
            
            [Mass0,K0,A0,F0,Q0,G0,H0,R0] = self.getDynamicFEMatrices(u0, t0);
            
            % Adjust the warnings.
            warnoffId = { 'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix', 'MATLAB:illConditionedMatrix'}; 
            for i = 1:length(warnoffId)    
              warnstat(i) = warning('query',warnoffId{i});
              warnoff(i) = warnstat(i);
              warnoff(i).state = 'off';
            end

            warning(warnoff);
            
            % Check for NaN propagation by passing t as NaN
            [Mass1,K1,A1,F1,Q1,G1,H1,R1] = self.getDynamicFEMatrices(u0, t1);
            
            % Check for NaN propagation in mass matrix by passing u as NaN
            u1 = NaN*ones(size(u0,1), 1);
            [Mass2,~,~,~,~,~,~,~] = self.getDynamicFEMatrices(u1, t0);
            
            warning(warnstat);
            
            if(any(isnan(Mass1(:))))
                self.dt=true;
            end
            
            if(~self.dt && any(isnan(Mass2(:))))
                self.du=true;
            end
            
            
            [self.ct,self.at,self.ft,...
                self.qt,self.gt,self.ht,self.rt] = self.checkNaNFEMatrix(K1,A1,F1,Q1,G1,H1,R1);
            
            self.vq = self.qt || self.qu;
            self.vg = self.gt || self.gu;
            self.vh = self.ht || self.hu;
            self.vr = self.rt || self.ru;
            self.vc = self.ct || self.cu;
            self.va = self.at || self.au;
            self.vf = self.ft || self.fu;
            %dependency of both 'm' and 'd' coefficients on either time or u.
            self.vd = self.dt || self.du;
            
            self.IsSpatialCoefficientsNonlinear = (self.vc || self.va || ...
                self.vf || self.vq || self.vg || self.vh || self.vr);
            
            % Update the system matrices.
            self.Mass = Mass0;
            self.K = K0;
            self.A = A0;
            self.F = F0;
            self.Q = Q0;
            self.G = G0;
            self.H = H0;
            self.R = R0;
            
            [self.Nu,self.Or]=pdenullorth(H0);
            if size(self.Or,2)==0
                ud0=zeros(size(K0,2),1);
            else
                ud0=full(self.Or*((H0*self.Or)\R0));
            end
            
            [self.totalNumEqns, self.numConstrainedEqns]=size(self.Nu);
            
            self.B=self.Nu;
            self.ud=ud0;
            self.dudt=0*ud0; %Assume zero initial derivatives w.r.t. time. 
            self.tspan=max(tlist)-min(tlist);
            
            
            % determine which nodes are within ignored regions
            % these will have undefined K
            ind = find(max(abs(K0),[],2) ~= 0);
            Nind = length(ind);
            self.Ni = sparse(ind,1:Nind,1,self.totalNumEqns,Nind);
            self.numIgnoredEqns = self.totalNumEqns-Nind;
            
            [Bind,~] = find(self.B);
            self.Bkeep = ismember(Bind,ind);
            self.B = self.B(:,self.Bkeep);
            
            % determine which nodes are not dynamic (infinitely stiff)
%             ind = find(max(abs(Mass0),[],2) ~= 0);



        end
        
        function self = checkQforUDependence(self,u0,tdummy)
            
            %checkFuncDepen determine whether coefficients are functions on u
            [~,~,~,Q0,~,~,~] = self.getStationaryFEMatrices(u0,tdummy);
            u1 = NaN*ones(size(u0,1), 1);
            [~,~,~,Q1,~,~,~] = self.getStationaryFEMatrices(u1,tdummy);
            
            % deal with outflow boundaries
            self.outflow_ind = any(isnan(Q0)')';
            Q1(self.outflow_ind,:) = 0;
            
            % check if coefficients are functions of u
            if any(isnan(Q1(:)))
                self.qu = true; 
            else
                self.qu = false; 
            end
            
        end
    end
    
    
    methods(Access=private)
        
        function [Mass,K,A,F,Q,G,H,R] = getDynamicFEMatrices(self, u, time)
            
            % Adjust the warnings.
            warnoffId = { 'MATLAB:singularMatrix', 'MATLAB:nearlySingularMatrix'}; 
            for i = 1:length(warnoffId)    
              warnstat(i) = warning('query',warnoffId{i});
              warnoff(i) = warnstat(i);
              warnoff(i).state = 'off';
            end
            
            warning(warnoff);
            
            [K,A,F,Q,G,H,R]  = self.getStationaryFEMatrices(u,time);
            
            if(self.nrp==2)
                [B,KB, FB, AB] = formGlobalSUPG2D(self.emptyPDEModel, self.p, self.t, self.coefstruct, self.coefstructb,u,time);
            elseif(self.nrp==3)
                [B,KB, FB, AB] = formGlobalSUPG3D(self.emptyPDEModel, self.p, self.t, self.coefstruct, self.coefstructb,u,time);
            end
            
            warning(warnstat);
            
            % deal with outflow boundaries
%             allQ_ind = sum(abs(Q),2) ~= 0;
            Q(self.outflow_ind,:) = 0;
            KB(self.outflow_ind,:) = 0;
            B(self.outflow_ind,:) = -B(self.outflow_ind,:);
            
            % add in petrov corrections
            K = K + B + KB; % b terms and Petrov corrections
            A = A + AB; % Petrov corrections
            F = F + FB; % Petrov corrections
            
            if(self.nrp==2)
                if(self.IsSecondOrderODE)
                    Mass = formGlobalM2D(self.emptyPDEModel, self.p, self.t, self.coefstruct,u,time,'m');
                else
                    Mass = formGlobalM2D(self.emptyPDEModel, self.p, self.t, self.coefstruct,u,time,'d');
                end
                
            elseif(self.nrp==3)
                if(self.IsSecondOrderODE)
                    Mass = formGlobalM3D(self.emptyPDEModel, self.p, self.t, self.coefstruct,u,time,'m');
                else
                    Mass = formGlobalM3D(self.emptyPDEModel, self.p, self.t, self.coefstruct,u,time,'d');
                end
            end
        end
    end
    
    methods(Access=protected, Static)
        function checkTlist(tlist)
            if length(tlist) < 2
                error(message('pde:pdeODEInfo:tlistTooShort'));
            end
            t0 = tlist(1);
            tfinal = tlist(end);
            if(t0 == tfinal)
                error(message('pde:pdeODEInfo:tlistEndpointsNotDistinct'));
            end
            tdir = sign(tfinal - t0);
            if any( tdir*diff(tlist) <= 0 )
                error(message('pde:pdeODEInfo:tlistUnordered'));
            end
        end
    end
    
end

