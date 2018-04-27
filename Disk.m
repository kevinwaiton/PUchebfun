classdef Disk < Domain
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius = 1;
        center = [0 0];
    end
    
    methods
        function obj = Disk(radius,center)
            obj.radius = radius;
            obj.center = center;
        end
        
        function ind = Interior(obj,pts)
            ind = sum((pts-obj.center).^2,2)<=obj.radius.^2;
        end
        
        
        function bound = Boundary(obj,N)
            
            THo = linspace(-pi,pi,N)';
            
            bound = obj.radius*[cos(THo)+obj.center(1) sin(THo)+obj.center(2)];
            
        end
        
        function plot(obj)
            
            Nxf = 200;
            Nyf = 200;
            
            x_f = linspace(obj.center(1)-obj.radius,obj.center(1)+obj.radius,Nxf)';
            y_f = linspace(obj.center(2)-obj.radius,obj.center(2)+obj.radius,Nyf)';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            ind = obj.Interior([Xf(:) Yf(:)]);
            
            P = [Xf(ind) Yf(ind)];
            
            TRI = delaunay(P);
            
            Z = zeros(length(P),1);
            
            trisurf(TRI,Xf(ind),Yf(ind),Z);
            view(0,90);
            shading interp;
        end
    end
end

