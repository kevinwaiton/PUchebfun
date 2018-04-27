classdef LSPatch3D < LSPatch
    % LSPatch2D PUFun class for representing n-d functions on non-square domains.
    %
    % This class represents a single tensor product polynomial, where the
    % domain of the Chebyshev polynomial bounds the domain of the function to
    % be approximated. The coefficients are found by solving a rank deficient
    % least square problem that minimizes the l_2 norm between a given function
    % f and the Chebyshev polynomial for a set of points inside the domain of
    % f.
    
    % LSPatch2D(varargin) constructs a tensor product approximation
    % representing a function, based on the options passed into with varargin;
    % that is PUFun('perf1',perf1,'pref2',pref2,..) is called. This
    % preferences that can be set are:
    %
    % The max lengths of the patches before sampling is to occur:
    % 'MaxLengths', [d_1 d_2 d_3]
    %
    % *The non square domain: 'InnerDomain', domain object
    %
    % *The domain used for the Chebyshev polynomial: 'domain', [a,b;c,d;e,f]
    %
    % *The zone (non overlapping part from partition) used: 'zone', [a,b;c,d;e,f]
    %
    % *The domain of the root of the tree: 'outerbox', [a,b;c,d;e,f]
    %
    % *An array of boolean indicies indicating if the approximation can be
    % split in a given dimension: 'canSplit', [bool_1,bool2]
    %
    % *The tolerance used for refinement: 'tol', 1e-b
    %
    % *The degree indices from the standard degrees in each dimension for non
    % square domains : 'degreeIndex', [ind_1,ind_2 ind_3].
    %
    % *The coarse degree to be used (if applicable)
    % : 'coarseDegreeIndex', [ind_1,ind_2 ind_3].
    %
    % *The degree indices from the standard degrees in each dimension for
    % square domains : 'ChebDegreeIndex', [ind_1,ind_2 ind_3].
    %
    % Here the degrees can be chosen from the set [3 5 9 17 33 65 129].
    % So if 'degreeIndex', [5 5 5], the max degree of any approximate will be
    % 33 in each direction.
    %
    % LSPatch2D(struct) will construct an approximation with a structure
    % struct. Here struct must contain the following fields:
    % in_domain : inner non square domain
    % outerbox : domain of the outerbox
    % zone : domain of the zone
    % domain : square domain of the polynomial
    % deg_in : indicies of degree for polynomials representing non square domains
    % cheb_deg_in : indicies of degree for polynomials representing square domains
    % cdeg_in : indicies of coarse degree
    % split_flag : boolean array indiciating to split in a dimension or not
    % max_lengths : obj.max_lengths: The max lengths of the patches before sampling is to occur
    % tol : tolerance used for refinement
    
    
    properties
        mid_values_err = inf %Store the evaluation at the Cheb points of the first kind
        MAT
    end
    
    properties (Access = protected)
        swap_deg_in
    end
    
    methods
        % function obj = LSPatch2D(in_domain,max_lengths,domain,zone,outerbox,deg_in,cheb_deg_in,split_flag,tol,cdeg_in)
        function obj = LSPatch3D(varargin)
            
            if length(varargin)==1
                varargin = varargin{:};
            end
            
            %Call superclass constructor
            obj = obj@LSPatch(varargin);
            
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
        
        % This method samples the polynomial with a LS square fit
        % of the domain.
        %
        %     Input: f, grid_opt,fast_opt
        %         f: objective function to be sampled
        %  grid_opt: option if takes a grid determine from vectors x y,
        %            i.e. [X,Y] = f({x y});
        %  fast_opt: if the function can be sampled everywhere (an
        %            example being a polynomial) then FFT is used.
        function max_val = sample(obj,f,grid_opt,~)
            
            if(nargin==2)
                grid_opt = false;
            end
            
            max_val = 0;
            
            
            x = chebpts(obj.degs(1)*2,obj.domain(1,:));
            y = chebpts(obj.degs(2)*2,obj.domain(2,:));
            z = chebpts(obj.degs(2)*2,obj.domain(3,:));
            
            [X,Y,Z] = ndgrid(x,y,z);
            
            XP = [X(:) Y(:) Z(:)];
            
            ind = obj.domain_in.Interior(XP);
            
            x_in = (1:obj.degs(1))';
            y_in = (1:obj.degs(2))';
            z_in = (1:obj.degs(3))';
            
            [X_in, Y_in, Z_in] = ndgrid(x_in,y_in,z_in);
            
            ind_c = sqrt(X_in.^2+Y_in.^2+Z_in.^2)<=max(obj.degs);
            
            Mx = chebtech.clenshaw(chebpts(obj.degs(1)*2),eye(obj.degs(1)));
            My = chebtech.clenshaw(chebpts(obj.degs(2)*2),eye(obj.degs(2)));
            Mz = chebtech.clenshaw(chebpts(obj.degs(3)*2),eye(obj.degs(3)));
            
            M = kron(Mz,kron(My,Mx));
            
            obj.MAT = M(ind,:);
            
            M = M(ind,:);
            
            if~ grid_opt
                F = f(X(ind),Y(ind),Z(ind));
            else
                F = f({x,y,z});
                F = F(ind);
            end
            
            max_val = max(abs(F));
            
            obj.coeffs = zeros(prod(obj.degs),1);
            
            warning('off','all');
            obj.coeffs(ind_c) = M\F;
            warning('on','all');
            
            
            obj.coeffs = reshape(obj.coeffs,obj.degs);
            
            E = obj.evalfGrid({x,y,z});
            E = E(ind);
            E = E(:) - F;
            
            %This is used to determin the point wise error
            obj.mid_values_err = max(abs(E(:)))./max(abs(F));
            
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
            
            if obj.mid_values_err>obj.tol
                
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
                %Chop(obj);
                Child = obj;
                Child.is_refined = true;
            end
            
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
            
            
            [new_zone_fit0,new_domain_fit0] = LSPatch3D.tightenLeaf(zone0,domain0,obj.outerbox,obj.domain_in);
            [new_zone_fit1,new_domain_fit1] = LSPatch3D.tightenLeaf(zone1,domain1,obj.outerbox,obj.domain_in);
            
            %zone0(split_dim,:) = new_zone_fit0(split_dim,:);
            %domain0(split_dim,:) = new_domain_fit0(split_dim,:);
            %
            %zone1(split_dim,:) = new_zone_fit1(split_dim,:);
            %domain1(split_dim,:) = new_domain_fit1(split_dim,:);
            
            zone0 = new_zone_fit0;
            domain0 = new_domain_fit0;
            
            zone1 = new_zone_fit1;
            domain1 = new_domain_fit1;
            %We first figure out if the the subdomains sit entirely in the domain itself.
            %In this case, we would just use a standard chebyshev
            %tensor product approximation.
            x1 = chebpts(16,domain0(1,:))';
            y1 = chebpts(16,domain0(2,:))';
            z1 = chebpts(16,domain0(3,:))';
            
            [X1,Y1,Z1] = ndgrid(x1,y1,z1);
            XP1 = [X1(:),Y1(:) Z1(:)];
            
            x2 = chebpts(16,domain1(1,:))';
            y2 = chebpts(16,domain1(2,:))';
            z2 = chebpts(16,domain1(3,:))';
            
            [X2,Y2,Z2] = ndgrid(x2,y2,z2);
            
            XP2 = [X2(:),Y2(:),Z2(:)];
            
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
                children{1} = LSPatch3D(struct0);
            end
            
            if all(obj.domain_in.Interior(XP2))
                %The square is in the domain. Set the child to a
                %standard Chebpatch
                children{2} = ChebPatch(struct1);
            else
                %The square is not in the domain. Set the child to a
                %least square patch
                children{2} = LSPatch3D(struct1);
            end
            
            x = chebpts(16,obj.domain(1,:))';
            y = chebpts(16,obj.domain(2,:))';
            z = chebpts(16,obj.domain(3,:))';
            
            [X,Y,Z] = ndgrid(x,y,z);
            
            XP = [X(:),Y(:),Z(:)];
            
            ind = obj.domain_in.Interior(XP);
            
            XP = XP(ind,:);
            
            ind11 = XP(:,split_dim)<=domain0(split_dim,2);
            ind22 = XP(:,split_dim)>=domain1(split_dim,1);
            
            
            if all(ind11)
                %The domain sits entirely in the first child
                Child = children{1};
                
            elseif all(ind22)
                
                %The domain sits entirely in the second child
                Child = children{2};
            else
                %Return the PUPatch with the new children
                Child = PUPatch(obj.domain,obj.zone,children,split_dim);
            end
            
            if set_vals
                for k=1:2
                    Child.children{k}.sample(obj.evalfGrid(Child.children{k}.leafGrids()));
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
        function [new_zone,new_domain] = splitleafGeom(zone,domain,outerbox,domain_in)
            
            obj.is_geometric_refined = true;
            
            x = linspace(domain(1,1),domain(1,2),50)';
            y = linspace(domain(2,1),domain(2,2),50)';
            z = linspace(domain(3,1),domain(3,2),50)';
            
            [X,Y,Z] = ndgrid(x,y,z);
            
            XP = [X(:) Y(:) Z(:)];
            
            ind = domain_in.Interior(XP);
            
            XP = XP(ind,:);
            
            new_domain = zeros(3,2);
            
            new_domain(:,1) = min(XP);
            
            new_domain(:,2) = max(XP);
            
            %pudge out zone a bit
            deltax = 0.5*Patch.overlap*diff(zone(1,:));

            deltay = 0.5*Patch.overlap*diff(zone(2,:));
            
            deltaz = 0.5*Patch.overlap*diff(zone(3,:));
            
            new_domain(1,1) = max(new_domain(1,1)-deltax,domain(1,1));
            new_domain(1,2) = min(new_domain(1,2)+deltax,domain(1,2));
            
            new_domain(2,1) = max(new_domain(2,1)-deltay,domain(2,1));
            new_domain(2,2) = min(new_domain(2,2)+deltay,domain(2,2));
            
            new_domain(3,1) = max(new_domain(3,1)-deltaz,domain(3,1));
            new_domain(3,2) = min(new_domain(3,2)+deltaz,domain(3,2));
           
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
            
            if abs(new_domain(3,1)-domain(3,1))>1e-10
                new_zone(3,1) = new_domain(3,1);
            end
            
            if abs(new_domain(3,2)-domain(3,2))>1e-10
                new_zone(3,2) = new_domain(3,2);
            end
        end
        end
end