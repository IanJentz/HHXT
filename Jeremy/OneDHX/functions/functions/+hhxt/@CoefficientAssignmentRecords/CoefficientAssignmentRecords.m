classdef (Sealed) CoefficientAssignmentRecords  < handle & matlab.mixin.internal.CompactDisplay
% CoefficientAssignmentRecords  Assignments of the equation coefficients
%    The CoefficientAssignmentRecords holds a record of all the assignments of 
%    equation coefficients across the geometric domain. This may be a single
%    assignment that represents the same field equations throughout the domain
%    or multiple assignments where the field equations are different across
%    various subdomains. This object cannot be created directly, it can only 
%    be created by the PDEModel class.
%
% CoefficientAssignmentRecords methods:
%    globalCoefficients - True if consistent across the geometric domain
%    findCoefficients   - Find coefficients assignment for a geometric region
%
% CoefficientAssignmentRecords properties:
%    CoefficientAssignments - Vector of coefficient assignments
%
% See also pde.PDEModel, pde.PDEModel/specifyCoefficients

% Copyright 2015-2017 The MathWorks, Inc.

properties
    
% CoefficientAssignments - Vector of coefficient assignments
%    A vector containing the coefficient assignments to the geometric domain
%    or subdomains. Each entry in the vector is a pde.CoefficientAssignment object. 
%    This object defines the coefficients of the PDE that should be imposed
%    on a designated domain or subdomain. Coefficients and host domains are 
%    defined using the specifyCoefficients method of the PDEModel.
    CoefficientAssignments;

end

methods
   % Method declaration
   tf = globalCoefficients(self, varargin)    
   ca = findCoefficients(self, varargin)   
   tf = coefficientsSpecified(self, varargin)  
end

methods %(Hidden=true, Access={?pde.EquationModel, ?pde.PDEResults})
    function obj=CoefficientAssignmentRecords(pdem) 
       obj.CoefficientAssignments = [];
       obj.ParentPdemodel = pdem;
    end  
    
    function tf = allCoefficientsNumeric(self)
        tf = true;
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients'));
        end
        if self.ParentPdemodel.IsTwoD
            qregiontype = 'face';
            nRegion = self.ParentPdemodel.Geometry.NumFaces;
        else
            qregiontype = 'cell';
            nRegion = self.ParentPdemodel.Geometry.NumCells;
        end
        
            for i = 1:nRegion
                ca = self.findCoefficients(qregiontype,i);
                if ~numericCoefficients(ca)
                    tf = false;
                    return
                end
            end
    end
    
    % Check that each subdomain has a coefficient assignment
    % also perform a basic sanity check on the coefficient
    % in the event that it was edited since assigned.
    % If both m and d are defined, then d must be a matrix.
    function performSolverPrechecks(self)
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients'));
        end
        if ~consistentCoefficientDefinition(self)
            error(message('pde:pdeCoefficientSpecification:inconsistentCoefPresence'));
        end               
        
        syssize = self.ParentPdemodel.PDESystemSize;
        if self.ParentPdemodel.IsTwoD
            qregiontype = 'face';
            nRegion = self.ParentPdemodel.Geometry.NumFaces;
        else
            qregiontype = 'cell';
            nRegion = self.ParentPdemodel.Geometry.NumCells;
        end
        nf = self.ParentPdemodel.Geometry.NumFaces;
        nc = self.ParentPdemodel.Geometry.NumCells;
        for i = 1:nRegion
            ca = self.findCoefficients(qregiontype,i);
            ca.performSolverPrecheck(syssize, nf, nc);
        end
    end
    
    function coefstruct = packCoefficients(self)
        [loc, state] = getSampleLocState(self);
        msh = self.ParentPdemodel.Mesh;
        ma = msh.MeshAssociation;
        numelems = size(msh.Elements,2);
        
        if self.globalCoefficients()
            coefstruct.ElementsInSubdomain = (1:numelems)';
            coefstruct.NumElementsInSubdomain = numelems;
            ca = self.CoefficientAssignments(end);
            thisstruct.b{1} = ca.b;
            thisstruct.d{1} = ca.d;
            thisstruct.c{1} = ca.c;
            thisstruct.a{1} = ca.a;
            coefstruct.Coefficients = thisstruct;
            coefstruct.IsComplex = hasComplexCoefficient(ca, loc, state);
        else
            if (size(msh.Nodes,1) == 2)
                qregiontype = 'face';
                nRegion = self.ParentPdemodel.Geometry.NumFaces;
                eidtofid = ma.RegionAssociativity; % Element ID to face ID.
                coefstruct.ElementsInSubdomain = zeros(numelems,1)';
                coefstruct.NumElementsInSubdomain = zeros(nRegion,1);
                complexcoef = false(4,1);
                endidx = 0;
                for i = 1:nRegion
                    felems = find(eidtofid==i);
                    numfelems = numel(felems);
                    startidx = endidx;
                    endidx = startidx+numfelems;
                    startidx = startidx+1;
                    coefstruct.NumElementsInSubdomain(i) = numfelems;
                    coefstruct.ElementsInSubdomain(startidx:endidx)=felems;
                    ca = self.findCoefficients(qregiontype,i);
                    
                    thisstruct.b{i} = ca.b;
                    thisstruct.d{i} = ca.d;
                    thisstruct.c{i} = ca.c;
                    thisstruct.a{i} = ca.a;
                    
%                     %test to see which coefficients evaluate
%                     region = struct('x',0,'y',0,'z',0,'subdomain',i);
%                     nan_n = nan*ones(self.ParentPdemodel.PDESystemSize,1);
%                     nan_d = nan*ones(size(msh.Nodes,1),1);
%                     state = struct('u',nan_n,'ux',nan_n,'uy',nan_n,'uz',nan_n,'time',0,...
%                         'kx',nan_d,'ky',nan_d,'kz',nan_d,...
%                         'uDx',nan_d,'uDy',nan_d,'uDz',nan_d,...
%                         'Rex',nan_d,'Rey',nan_d,'Rez',nan_d,...
%                         'uDDELTx',nan_d,'uDDELTy',nan_d,'uDDELTz',nan_d);
%                     
%                     if isa(thisstruct.b{i},'function_handle') 
%                         if ~any(isnan(thisstruct.b{i}(region,state))); thisstruct.b{i} = 0; end
%                     end
%                      if isa(thisstruct.d{i},'function_handle') 
%                         if ~any(isnan(thisstruct.d{i}(region,state))); thisstruct.d{i} = 0; end
%                      end
%                      if isa(thisstruct.c{i},'function_handle') 
%                         if ~any(isnan(thisstruct.c{i}(region,state))); thisstruct.c{i} = 0; end
%                      end
%                      if isa(thisstruct.a{i},'function_handle') 
%                         if ~any(isnan(thisstruct.a{i}(region,state))); thisstruct.a{i} = 0; end
%                      end

                    coefstruct.Coefficients = thisstruct;
                    complexcoef = complexcoef | hasComplexCoefficient(ca, loc, state);
                end
            else
                qregiontype = 'cell';
                nRegion = self.ParentPdemodel.Geometry.NumCells;
                eIDtocellID = ma.SubdomainAssociativity;
                coefstruct.ElementsInSubdomain = zeros(numelems,1)';
                coefstruct.NumElementsInSubdomain = zeros(nRegion,1);
                complexcoef = false(5,1);
                endidx = 0;
                for i = 1:nRegion
                    celems = eIDtocellID(2,i)-eIDtocellID(1,i)+1;
                    numcelems = celems;
                    startidx = endidx;
                    endidx = startidx+numcelems;
                    startidx = startidx+1;
                    coefstruct.NumElementsInSubdomain(i) = numcelems;
                    coefstruct.ElementsInSubdomain(startidx:endidx)=eIDtocellID(1,i):eIDtocellID(2,i);
                    ca = self.findCoefficients(qregiontype,i);
                    
                    thisstruct.b{i} = ca.b;
                    thisstruct.d{i} = ca.d;
                    thisstruct.c{i} = ca.c;
                    thisstruct.a{i} = ca.a;
                    
%                     %test to see which coefficients evaluate
%                     region = struct('x',0,'y',0,'z',0,'subdomain',i);
%                     nan_n = nan*ones(self.ParentPdemodel.PDESystemSize,1);
%                     nan_d = nan*ones(size(msh.Nodes,1),1);
%                     state = struct('u',nan_n,'ux',nan_n,'uy',nan_n,'uz',nan_n,'time',0,...
%                         'kx',nan_d,'ky',nan_d,'kz',nan_d,...
%                         'uDx',nan_d,'uDy',nan_d,'uDz',nan_d,...
%                         'Rex',nan_d,'Rey',nan_d,'Rez',nan_d,...
%                         'uDDELTx',nan_d,'uDDELTy',nan_d,'uDDELTz',nan_d);
%                     
%                     if isa(thisstruct.b{i},'function_handle') 
%                         if ~any(isnan(thisstruct.b{i}(region,state))); thisstruct.b{i} = 0; end
%                     end
%                      if isa(thisstruct.d{i},'function_handle') 
%                         if ~any(isnan(thisstruct.d{i}(region,state))); thisstruct.d{i} = 0; end
%                      end
%                      if isa(thisstruct.c{i},'function_handle') 
%                         if ~any(isnan(thisstruct.c{i}(region,state))); thisstruct.c{i} = 0; end
%                      end
%                      if isa(thisstruct.a{i},'function_handle') 
%                         if ~any(isnan(thisstruct.a{i}(region,state))); thisstruct.a{i} = 0; end
%                      end
                     
                    coefstruct.Coefficients = thisstruct;
                    complexcoef = complexcoef | hasComplexCoefficient(ca, loc, state);
                    
                end
            end
            coefstruct.IsComplex = complexcoef;
        end
    end
end
    
methods(Hidden=true, Access={?pde.EquationModel,?pde.InitialConditionsRecords})    
    
    function tf = timeDependent(self)
        tf = (mDefined(self) || dDefined(self));
    end
end
  
methods(Hidden=true, Access=private)
    % Returns true if each coef has consistent definition across all
    % subdomains. For example, m coefficient defined in one region is
    % also defined in all other regions.
    function tf = consistentCoefficientDefinition(self)
        tf = globalCoefficients(self);
        if tf
            return
        end
        tf = true;
        if self.ParentPdemodel.IsTwoD
            qregiontype = 'face';
            nRegion = self.ParentPdemodel.Geometry.NumFaces;
        else
            qregiontype = 'cell';
            nRegion = self.ParentPdemodel.Geometry.NumCells;
        end

        thiscoef = self.findCoefficients(qregiontype, 1);
        for i = 2:nRegion
            othercoef = self.findCoefficients(qregiontype, i);           
            if ~thiscoef.coefficientsMatch(othercoef)
                tf = false;
                return;
            end            
        end                   
    end           
    
    function tf = hasComplexCoefficients(self)
        
        [location, state] = getSampleLocState(self);                       
        tf = false;
        if ~self.coefficientsSpecified()
            error(message('pde:pdeCoefficientSpecification:subdomainsWithoutCoefficients')); 
        end     
        if globalCoefficients(self)
            % Just need to check one.
            thiscoef = self.CoefficientAssignments(end);
            tf = any(thiscoef.hasComplexCoefficient(location, state));
            return
        end
        
        
        if self.ParentPdemodel.IsTwoD
            qregiontype = 'face';
            nRegion = self.ParentPdemodel.Geometry.NumFaces;
        else
            qregiontype = 'cell';
            nRegion = self.ParentPdemodel.Geometry.NumCells;
        end
        
        for i = 1:nRegion
            thiscoef = self.findCoefficients(qregiontype, i);           
            if any(thiscoef.hasComplexCoefficient(location, state))
                tf = true;
                return;
            end            
        end                   
    end 
    
    function [location, state] = getSampleLocState(self)
        msh = self.ParentPdemodel.Mesh;
        nodecoords =msh.Nodes(:,1);
        location.x = nodecoords(1);
        location.y = nodecoords(2);
        if numel(nodecoords) == 3
            location.z = nodecoords(3);
        else
            location.z = 0;
        end
        location.subdomain=1;         
        N = self.ParentPdemodel.PDESystemSize;
        state.u(1:N,1)=0;   % state.u(1:NSystemSize, NumLocations)
        state.ux(1:N,1)=0;
        state.uy(1:N,1)=0;
        state.uz(1:N,1)=0;
        state.time=0;        
         str = round((N-1)/2);
        state.kx = zeros(str, 1);
        state.ky = zeros(str, 1);
        state.kz = zeros(str, 1);
        state.uDx = zeros(str, 1);
        state.uDy = zeros(str, 1);
        state.uDz = zeros(str, 1);
        state.Rex = zeros(str, 1);
        state.Rey = zeros(str, 1);
        state.Rez = zeros(str, 1);
        state.uDDELTx = zeros(str, 1);
        state.uDDELTy = zeros(str, 1);
        state.uDDELTz = zeros(str, 1);
    end
    
end

methods(Hidden=true, Access={?pde.CoefficientAssignment})
    function delistCoefficientAssignment(self, coeftodelete)
        numcoef = numel(self.CoefficientAssignments);
        for i = 1:numcoef
            thiscoef = self.CoefficientAssignments(i);
            if thiscoef == coeftodelete
                self.CoefficientAssignments(i) = [];
                break
            end
        end  
        numcoef = numel(self.CoefficientAssignments);
        if numcoef == 0 && isvalid(self.ParentPdemodel)
           self.ParentPdemodel.delistCoefficientAssignments(); 
        end
    end 
end

methods(Static, Hidden=true)    
    % preDelete - Preserved for backward compatibility. In R2017b and prior
    % releases, the coefficient assignment has an ObjectBeingDestroyed 
    % callback that calls this function.
    function preDelete(self,~)
        if isvalid(self.ParentPdemodel)
            self.ParentPdemodel.delistCoefficientAssignments();
        end
    end             
end

methods(Hidden=true)
    function delete(self)
       numcoef = numel(self.CoefficientAssignments);
        for i = 1:numcoef       
            if isvalid(self.CoefficientAssignments(i))           
                delete(self.CoefficientAssignments(i));
            end
        end            
        if isvalid(self.ParentPdemodel)
            self.ParentPdemodel.delistCoefficientAssignments();
        end      
    end     
end

properties (Hidden = true, SetAccess='private')
    ParentPdemodel;
end  
  
end
