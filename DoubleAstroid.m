classdef DoubleAstroid < Domain
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = DoubleAstroid()
        end
        
        function ind = Interior(obj,pts)
            
            [THo,R] = cart2pol(pts(:,1),pts(:,2));
            
            TH = mod(THo-pi+pi/8,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            ind1 = R<=(cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2);
           
            TH = mod(THo-pi-pi/8,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            ind2 = R<=(cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2);
            
            ind = ind1 | ind2;
            
        end
        
        function bound = Boundary(obj,N)
            
            THo = linspace(-pi,pi,N)';
            
            TH = mod(THo-pi+pi/8,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            bound_r_1 = (cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2);
            
            TH = mod(THo-pi-pi/8,2*pi)-pi;
            
            TH = mod(TH,pi/2);
            
            bound_r_2 = [(cos(TH).^(2/3)+sin(TH).^(2/3)).^(-3/2)];
            
            [x,y] = pol2cart(THo,max(bound_r_1,bound_r_2));
            
            bound = [x(:) y(:)];
            
        end
        
        function plot(obj)
            
            Nxf = 200;
            Nyf = 200;
            
            x_f = linspace(-1,1,Nxf)';
            y_f = linspace(-1,1,Nyf)';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            B = obj.Boundary(800);
            
            ind = obj.Interior([Xf(:) Yf(:)]);
            
            P = [Xf(ind) Yf(ind)];
            
            TRI = delaunayTriangulation([B;P],[(1:length(B)-1)' (2:length(B))'; length(B) 1]);
            
            TF = isInterior(TRI);
            
            Z = zeros(length(P)+length(B),1);
            
            trisurf(TRI.ConnectivityList(TF,:),TRI.Points(:,1),TRI.Points(:,2),Z);
            view(0,90);
            shading interp;
        end
    end
end

