classdef (Abstract) Domain
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        boundary_points % boundary points. Can be empty
    end
    
    methods (Abstract)
        ind = Interior(obj,pts) %Give the indices of the pts in the domain
        bound = Boundary(obj,N) %Gives the boundary
    end
   
end

