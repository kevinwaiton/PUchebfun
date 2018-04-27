classdef LeafPatch<Patch
    % This is the abstract class for a leaf object. This is used with the
    % PUPatch object.
    
    % LeafPatch(var_struct) serves as the base constructor for a leafPatch
    % object, where var_struct is a structure with fields:
    %
    % *var_struct.domain, the domain used for the Chebyshev polynomial.
    %
    % *var_struct.zone, the zone (non overlapping part from partition). If the
    % zone is undefined, the zone is set to obj.domain.
    %
    % *var_struct.outerbox, the domain of the function. If the outerbox is
    % undefined, the outerbox is set to obj.domain.
    properties
        index = [];
        chebweights = [];
        cheb_bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
        bump
        values
        coeffs
        is_interp
        deg_in %index for the standard degrees
        degs %array of degrees along the dimensions
    end
    
    properties (Constant)
        standard_variables = load('cheb_points_matrices.mat');
        standard_degs = [3 5 9 17 33 65 129];
        invf = @(x,dom) 2/diff(dom)*x-sum(dom)/diff(dom); %takes points from a domain to [-1 1]
        forf = @(x,dom) 0.5*diff(dom)*x+0.5*sum(dom); %takes points from [-1 1] to a domain
    end
    
    methods
        function obj = LeafPatch(var_struct)
            
            if isfield(var_struct, 'domain')
                obj.domain = var_struct.domain;
            else
                error('A domain needs to be specified.');
            end
            
            if isfield(var_struct, 'outerbox')
                obj.outerbox = var_struct.outerbox;
            else
                obj.outerbox = obj.domain;
            end
            
            if isfield(var_struct, 'zone')
                obj.zone = var_struct.zone;
            else
                obj.zone = obj.domain;
            end
            
            
            obj.is_geometric_refined = true; %square is always geometrically refined
            [obj.dim,~] = size(obj.domain);
            obj.is_leaf = true;
            
            obj.bump = cell(3,1);
            
            for k=1:obj.dim
                % if isequal(obj.domain(k,:),obj.outerbox(k,:))
                %     w = @(x) ones(size(x));
                if obj.domain(k,1) == obj.outerbox(k,1)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))+1)/2);
                elseif obj.domain(k,2) == obj.outerbox(k,2)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))-1)/2);
                else
                    w = @(x) obj.cheb_bump(obj.invf(x,obj.domain(k,:)));
                end
                obj.bump{k} = w;
            end
        end
        
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: array of length(X) of the patch's weight evaluated at X.
        function ef = evalfBump(obj,X)
            
            ef = ones(size(X,1),1);
            for k=1:obj.dim
                ef = ef.*obj.bump{k}(X(:,k));
            end
            
        end
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: array of dim(X) of the patch's weight evaluated at X.
        function ef = evalfGridBump(obj,X)
            
            W = cell(3,1);
            
            for i=1:obj.dim
                W{i} = obj.bump{i}(X{i});
            end
            
            if obj.dim==1
                ef = W{1};
            elseif obj.dim==2
                ef = W{1}*W{2}.';
            else
                ef = reshape(W{3},1,1,length(W{3})).*(W{2}'.*W{1});
            end
            
        end
        
        % Plots the overlaping domain of the tree. Works only
        % in 2D and 3D.
        %
        %     Input:
        %     color: color of marker for the center of the domain.
        %
        function plotdomain(obj,color)
            
            if nargin==1
                color = 'black';
            end
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:))];
                rectangle('position',[obj.domain(:,1)' lengths'],'LineWidth',2,'EdgeColor',color);
                plot(mean(obj.domain(1,:)),mean(obj.domain(2,:)),'.','MarkerSize',10,'Color',color);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:));diff(obj.domain(3,:))];
                center = sum(obj.domain,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                plot3(mean(obj.domain(1,:)),mean(obj.domain(2,:)),mean(obj.domain(3,:)),'.','MarkerSize',10,'Color','black');
                hold off;
            end
        end
        
        
        % Plots the zones of the tree. Works only
        % in 2D and 3D.
        %
        function plotzone(obj)
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:))];
                rectangle('position',[obj.zone(:,1)' lengths'],'LineWidth',2);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:));diff(obj.zone(3,:))];
                center = sum(obj.zone,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                hold off;
            end
        end
        
        % Returns the points of the function
        function pts = points(obj)
            
            C = cell(obj.dim,1);
            
            for i=1:obj.dim
                C{i} = obj.standard_variables.chebpoints{obj.deg_in(i)};
                C{i} = obj.forf(C{i},obj.domain(i,:));
            end
            [out{1:obj.dim}] = ndgrid(C{:});
            
            pts = zeros(numel(out{1}),obj.dim);
            
            for i=1:obj.dim
                pts(:,i) = out{i}(:);
            end
            
        end
        
        % Returns a cell array of the grids of the domain
        function grid = leafGrids(obj)
            grid = cell(1,obj.dim);
            
            for i=1:obj.dim
                grid{i} = obj.standard_variables.chebpoints{obj.deg_in(i)};
                grid{i} = obj.forf(grid{i},obj.domain(i,:));
            end
        end
        
        % TODO. Figure out what to do here!
        function ln=length(obj)
            ln = prod(obj.degs);
        end
        
    end
    methods (Abstract)
        %This method will split the child, creating a new PUPatch. If the
        %obj does not need to split, the method returns obj.
        Child = splitleaf(obj);
    end
    
end
