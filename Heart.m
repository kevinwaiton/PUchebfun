classdef Heart
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mult = 5;
    end
    
    methods
        function obj = Heart()
        end
        
        function ind = Interior(obj,pts)
            [theta,rho] = cart2pol(pts(:,1),pts(:,2));
            ind = rho<= obj.mult*sin(theta)-sin(obj.mult*theta);
        end
        
        function plot(obj)
            
            Nxf = 66;
            Nyf = 66;
            
            x_f = chebpts(Nxf,[-3.1,3.1])';
            y_f = chebpts(Nyf,[0,4.6])';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            XPf = [Xf(:) Yf(:)];
            
            grid_sq_ind = obj.Interior(XPf);
            
            %fine grid in domain
            XPf = XPf(grid_sq_ind,:);
            
            scatter(XPf(:,1),XPf(:,2),'red');
        end
    end
    
end

