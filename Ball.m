classdef Ball < Domain
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius = 1;
        center = [0 0 0];
    end
    
    methods
        function obj = Ball(radius,center)
            obj.radius = radius;
            obj.center = center;
        end
        
        function ind = Interior(obj,pts)
            ind = sum((pts-obj.center).^2,2)<=obj.radius.^2;
        end
        
        function bound = Boundary(obj,N)
            bound = [];
        end
        
        function plot(obj)
            Nxf = 30;
            Nyf = 30;
            Nzf = 30;
            
            x_f = linspace(obj.center(1)-obj.radius,obj.center(1)+obj.radius,Nxf)';
            y_f = linspace(obj.center(2)-obj.radius,obj.center(2)+obj.radius,Nyf)';
            z_f = linspace(obj.center(3)-obj.radius,obj.center(3)+obj.radius,Nzf)';
            
            [Xf,Yf,Zf] = ndgrid(x_f,y_f,z_f);

            ind = obj.Interior([Xf(:) Yf(:) Zf(:)]);
            
            P = [Xf(ind) Yf(ind) Zf(ind)];
            
            TRI = delaunay(P);
            
            C = zeros(length(P),1);
             
            trisurf(TRI,Xf(ind),Yf(ind),Zf(ind),C);
            
        end
    end
end

