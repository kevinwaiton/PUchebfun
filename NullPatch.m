classdef NullPatch<Patch
% ChebPatch PUFun class for representing empty patches. While inelegant,
% it gives a simple way of intrucing null values for patches without
% breaking the rest of the code in PUPatch.
    properties
        values = 0;
        coeffs = 0;
        linOp
        ClinOp
    end
    
    methods
        % Construct for the ChebPatch
        %
        %      Input:
        %       zone: (dim x 2) array indiciating array for the zone.
        %     domain: (dim x 2) array indiciating array for the domain.
        %   outerbox: (dim x 2) array indiciating array for the outerbox.
        %    degs_in: (dim x 1) integer array indicating the degree in each
        %                       dimension. Here the standard degrees are
        %                       [3 5 9 17 33 65 129].
        % split_flag: (dim x 1) boolean array indicating if the patch can
        %                       be split in any given dimension.
        %function obj = ChebPatch(domain,zone,outerbox,deg_in,split_flag,tol,cdeg_in)
        function obj = NullPatch(var_struct)
            
            obj.is_leaf = true;
            obj.outerbox = var_struct.outerbox;
            obj.zone = var_struct.zone;
            obj.domain = var_struct.domain;
            dim = size(obj.domain,1);
            
            obj.split_flag = false(dim,1);
            obj.is_refined = true;
            obj.is_geometric_refined = true;
            obj.is_null = true;
            
        end
        
        % Returns the length of the object
        function ln=length(obj)
            ln = 0;
        end
        
        % Returns the values of the object
        function vals = Getvalues(obj)
            vals = 0;
        end
        
        
        % Returns the points of the function
        function pts = points(obj)
            pts = [];
        end
        
        % Returns a cell array of the grids of the domain
        function grid = leafGrids(obj)
            grid = [];
        end
        
        function sample(obj,f)
            
        end
        % Evaluates the approximant and its derivatives.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: length(X) array containing the interpolated
        function ef = evalf(obj,X,G)
            ef = 0;
        end
        
        % Evaluates the approximant on a grid.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: matrix of dim(X) containing the interpolated values
        function ef = evalfGrid(obj,X,G)
            ef = 0;
        end
        
        function plotdomain(obj)
        end
        
        
        end
        
        
    
end