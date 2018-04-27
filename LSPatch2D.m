classdef LSPatch2D < LSPatch
    % LSPatch2D PUFun class for representing n-d functions on non-square domains.
    %
    % This class represents a single tensor product polynomial, where the
    % domain of the Chebyshev polynomial bounds the domain of the function to
    % be approximated. The coefficients are found by solving a rank deficient
    % least square problem that minimizes the l_2 norm between a given function
    % f and the Chebyshev polynomial for a set of points inside the domain of
    % f.
    
    % LSPatch2D(var_struct) constructs a tensor product approximation
    % representing a function, based on the options passed into with
    % structure var_struct. The options are the same as LSPatch.m.
    
    properties
        mid_values_err = inf %Store the evaluation at the Cheb points of the first kind
        mid_values_err1 = inf;
        mult;
        
    end
    
    properties (Access = protected)
        swap_deg_in
    end
    
    
    
    methods
        % function obj = LSPatch2D(in_domain,max_lengths,domain,zone,outerbox,deg_in,cheb_deg_in,split_flag,tol,cdeg_in)
        function obj = LSPatch2D(var_struct)
            
            obj = obj@LSPatch(var_struct);
            
            obj.is_geometric_refined = true;
        end
        
        %Returns structure of parameters
        function p_struct = params(obj)
            
            p_struct.outerbox = obj.outerbox;
            p_struct.zone = obj.zone;
            p_struct.domain = obj.domain;
            p_struct.deg_in = obj.deg_in;
            p_struct.cheb_deg_in = obj.cheb_deg_in;
            p_struct.domain_in = obj.domain_in;
            p_struct.split_flag = obj.split_flag;
            p_struct.max_lengths = obj.max_lengths;
            p_struct.tol = obj.tol;
            p_struct.cdeg_in = obj.cdeg_in;
        end
    
        
        % This method determines if the leaf needs to be simplified,
        % and splits if necessary 
        %     Input:
        %       Max: Global max value
        %
        %    Output:
        %     Child: Returns the leaf if no splitting is needed, otherwise
        %            returns a PUPatch with the split leaf.
        function Child = splitleaf(obj,Max,~)
            
            obj.GlobalMax = Max;
            
            if obj.mid_values_err1>obj.tol
                
                Child = obj;
                %Go through and split in each unresolved direction
                for k=1:obj.dim
                    if Child.is_leaf
                        Child = split(Child,k);
                    else
                        Child.split(k);
                    end
                end
                
            else
                Chop(obj);
                Child = obj;
                Child.is_refined = true;
            end
            
        end
        
        
        
        % This method samples the polynomial with a LS square fit
        % of the domain.
        %
        %     Input: f, grid_opt,fast_opt
        %         f: objective function to be sampled
        %  grid_opt: option if takes a grid determine from vectors x y,
        %            i.e. [X,Y] = f({x y});
        %  fast_opt: if the function can be sampled everywhere (an
        %            example being a polynomial) then FFT is used.
        function max_val = sample(obj,f,grid_opt,fast_opt)
            
            
            
            if nargin==2
                grid_opt = false;
                fast_opt = false;
            elseif nargin==2
                fast_opt = false;
            end
            
            max_val = 0;
            
                
                for i=2:1:10
                    
                    mult = i;
                    
                    x = chebpts(obj.degs(1)*mult,obj.domain(1,:));
                    y = chebpts(obj.degs(2)*mult,obj.domain(2,:));
                    

                    
                    [X,Y] = ndgrid(x,y);
                    
            
                    
                    XP = [X(:) Y(:)];
                    

                    
                    ind = obj.domain_in.Interior(XP);

                    
                    if sum(ind)/prod(obj.degs)>=4
                        break;
                    end
                end
                
                x1 = chebpts(2*obj.degs(1),obj.domain(1,:),1);
                y1 = chebpts(2*obj.degs(2),obj.domain(2,:),1);
                
                [X1,Y1] = ndgrid(x1,y1);
                
                XP1 = [X1(:) Y1(:)];
                
                
                ind1 = obj.domain_in.Interior(XP1);
                
                obj.mult = mult;
                
                
                
                if ~grid_opt
                    F = f(X(ind),Y(ind));
                else
                    F = f({x,y});
                    F = F(ind);
                end
                
                if ~fast_opt
                    
                    Mx = chebtech.clenshaw(chebpts(obj.degs(1)*mult),eye(obj.degs(1)));
                    My = chebtech.clenshaw(chebpts(obj.degs(2)*mult),eye(obj.degs(2)));
                    
                    M = kron(My,Mx);
                    M = M(ind,:);
                    obj.coeffs = zeros(prod(obj.degs),1);
                    warning('off','all');
                    obj.coeffs = M\F;
                    warning('on','all');
                    
                    obj.coeffs = reshape(obj.coeffs,obj.degs);
                    
                else
                    if ~grid_opt
                        points = num2cell(obj.points(),1);
                        obj.coeffs = chebfun2.vals2coeffs(reshape(f(points{:}),obj.degs));
                    else
                        obj.coeffs = chebfun2.vals2coeffs(f(obj.leafGrids()));
                    end
                end
                
                
                if~ grid_opt
                    F1 = f(X1,Y1);
                else
                    F1 = f({x1,y1});
                end
                
                max_val = max(abs(F));
                
                obj.LocalMax = max_val;
                
                E = obj.evalfGrid({x,y});
                E = E(ind);
                E = E(:) - F;
                
                E1 = obj.evalfGrid({x1,y1});
                E1 = E1 - F1;
                E1 = E1(ind1);
                
                %This is used to determin the point wise error
                obj.mid_values_err = max(abs(E(:)))./max(abs(F)*sqrt(sum(ind)));
                obj.mid_values_err1 = max(abs(E1(:)))./max(abs(F1(ind1))*sqrt(sum(ind1)));
                
        end
        
        % The method determines will split a child into along
        % a dimension.
        %
        %     Input:
        %   overlap: overlap intended to be used for the splitting
        %
        %    Output:
        %     Child: the PUPatch with the two new children.
        function Child = split(obj,split_dim,set_vals)
            
            if nargin == 2
                set_vals = false;
            end
            
            children = cell(1,2);
            
            %The width of the overlap
            delta = 0.5*obj.overlap*diff(obj.zone(split_dim,:));
            
            zone0 = obj.zone;
            zone1 = obj.zone;
            
            domain0 = obj.domain;
            domain1 = obj.domain;
            
            m = sum(obj.zone(split_dim,:))/2;
            
            zone0(split_dim,:) = [obj.zone(split_dim,1) m];
            zone1(split_dim,:) = [m obj.zone(split_dim,2)];
            
            domain0(split_dim,:) = [max(obj.outerbox(split_dim,1),obj.zone(split_dim,1)-delta) m+delta];
            domain1(split_dim,:) = [m-delta,min(obj.outerbox(split_dim,2),obj.zone(split_dim,2)+delta)];
            
           [new_zone_fit0,new_domain_fit0] = LSPatch2D.tightenLeaf(zone0,domain0,obj.outerbox,obj.domain_in);
           [new_zone_fit1,new_domain_fit1] = LSPatch2D.tightenLeaf(zone1,domain1,obj.outerbox,obj.domain_in);
            
            
            zone0 = new_zone_fit0;
            domain0 = new_domain_fit0;
            
            zone1 = new_zone_fit1;
            domain1 = new_domain_fit1;
            
            %We first figure out if the the subdomains sit entirely in the domain itself.
            %In this case, we would just use a standard chebyshev
            %tensor product approximation.
            x1 = chebpts(16,domain0(1,:))';
            y1 = chebpts(16,domain0(2,:))';
            [X1,Y1] = ndgrid(x1,y1);
            XP1 = [X1(:),Y1(:)];
            
            x2 = chebpts(16,domain1(1,:))';
            y2 = chebpts(16,domain1(2,:))';
            [X2,Y2] = ndgrid(x2,y2);
            XP2 = [X2(:),Y2(:)];
            
            struct0 = obj.params;
            struct0.domain = domain0; struct0.zone = zone0;
            
            struct1 = obj.params;
            struct1.domain = domain1; struct1.zone = zone1;
            
            if all(obj.domain_in.Interior(XP1))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{1} = ChebPatch(struct0);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{1} = LSPatch2D(struct0);
                children{1}.is_geometric_refined = false;
            end
            
            if all(obj.domain_in.Interior(XP2))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{2} = ChebPatch(struct1);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{2} = LSPatch2D(struct1);
                children{2}.is_geometric_refined = false;
            end
            
            x = chebpts(16,obj.domain(1,:))';
            y = chebpts(16,obj.domain(2,:))';
            
            [X,Y] = ndgrid(x,y);
            
            XP = [X(:),Y(:)];
            
            ind = obj.domain_in.Interior(XP);
            
            XP = XP(ind,:);
            
            ind11 = XP(:,split_dim)<=domain0(split_dim,2);
            ind22 = XP(:,split_dim)>=domain1(split_dim,1);
            
            
            if all(ind11)
                %The domain sits entirely in the first child
                %Child = children{1};
                children{2} = NullPatch(struct1);
                Child = PUPatch(obj.domain,obj.zone,children,split_dim);
            elseif all(ind22)
                
                %The domain sits entirely in the second child
                %Child = children{2};
                children{1} = NullPatch(struct0);
                Child = PUPatch(obj.domain,obj.zone,children,split_dim);
            else
                %Return the PUPatch with the new children
                Child = PUPatch(obj.domain,obj.zone,children,split_dim);
                
                Child.is_geometric_refined = false;
            end
            
            
            if set_vals
                for k=1:2
                    Child.children{k}.sample(obj.evalfGrid(Child.children{k}.leafGrids()));
                end
            end
        end
        
        %This method attempts to simplify the polynomial
        function Chop(obj)
            
            loc_tol = obj.tol^(7/8);
            
            for k=1:obj.dim
                
                if obj.split_flag(k)
                    
                    colChebtech = chebfun3t.unfold(obj.coeffs, k);
                    colChebtech = sum(abs(colChebtech),2);
                    fCol = chebtech2({[],colChebtech});
                    hscale = diff(obj.domain(k,:));
                    
                    tol = loc_tol*max(obj.GlobalMax/obj.LocalMax,hscale);
                    cutoff = length(simplify(fCol, tol))+1;
                    
                    if cutoff<obj.degs(k)
                        obj.split_flag(k) = false;
                        j = find(cutoff<=obj.standard_degs);
                        
                        if j<obj.deg_in(k)
                            obj.deg_in(k) = j;
                            obj.deg(k) = obj.standard_degs(j);
                            
                            if k==1
                                obj.coeffs = obj.coeffs(1:obj.deg(k),:);
                            else
                                obj.coeffs = obj.coeffs(:,1:obj.deg(k));
                            end
                        end
                    end
                end
                
            end
            
        end
        
    end
    
    methods (Static)
        
        % This method fits a given square zone and domain to the inner
        % domain.
        %
        %     Input:
        %      zone,domain,outerbox: square zone, domain, and outerbox
        %                 domain_in: non-square inner domain to fit too.
        %    Output:
        %       new_zone,new_domain: fitted squared zone and domains.
        function [new_zone,new_domain] = tightenLeaf(zone,domain,outerbox,domain_in)
            
            obj.is_geometric_refined = true;
            
            x = linspace(domain(1,1),domain(1,2),240)';
            y = linspace(domain(2,1),domain(2,2),240)';
            
            [X,Y] = ndgrid(x,y);
            
            XP = [X(:) Y(:)];
            
            ind = domain_in.Interior(XP);
            
            XP = XP(ind,:);
            
            new_domain = zeros(2,2);
            
            new_domain(:,1) = min(XP);
            
            new_domain(:,2) = max(XP);
            
            %pudge out zone a bit
            deltax = 0.5*Patch.overlap*diff(zone(1,:));
            %The width of the overlap
            deltay = 0.5*Patch.overlap*diff(zone(2,:));
            
            new_domain(1,1) = max(new_domain(1,1)-deltax,domain(1,1));
            new_domain(1,2) = min(new_domain(1,2)+deltax,domain(1,2));
            
            new_domain(2,1) = max(new_domain(2,1)-deltay,domain(2,1));
            new_domain(2,2) = min(new_domain(2,2)+deltay,domain(2,2));
           
            new_zone = zone;
            
            if abs(new_domain(1,1)-domain(1,1))>1e-10
                new_zone(1,1) = new_domain(1,1);
            end
            
            if abs(new_domain(1,2)-domain(1,2))>1e-10
                new_zone(1,2) = new_domain(1,2);
            end
            
            if abs(new_domain(2,1)-domain(2,1))>1e-10
                new_zone(2,1) = new_domain(2,1);
            end
            
            if abs(new_domain(2,2)-domain(2,2))>1e-10
                new_zone(2,2) = new_domain(2,2);
            end
            
        end
    end
end