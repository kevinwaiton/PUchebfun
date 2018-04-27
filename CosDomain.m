classdef CosDomain
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        TimeLength;
        period;
        MinSpaceLength;
        MaxSpaceLength;
    end
    
    methods
        
        function obj = CosDomain(TimeLength,period,MinSpaceLength,MaxSpaceLength)
            obj.TimeLength = TimeLength;
            obj.period = period;
            obj.MinSpaceLength = MinSpaceLength;
            obj.MaxSpaceLength = MaxSpaceLength;
        end
        
        function ind = Interior(obj,pts)
            a = (obj.MaxSpaceLength+obj.MinSpaceLength)/2;
            b = (obj.MaxSpaceLength-obj.MinSpaceLength)/2;
            ind = pts(:,1)>=0 & pts(:,1) <= a+b*cos(pi/obj.period*pts(:,2));
        end
        
        function plot(obj)
            
            Nxf = 128;
            Nyf = 128;
            
            x_f = chebpts(Nxf,[0,obj.MaxSpaceLength])';
            y_f = chebpts(Nyf,[0,obj.TimeLength])';
            
            [Xf,Yf] = ndgrid(x_f,y_f);
            
            XPf = [Xf(:) Yf(:)];
            
            grid_sq_ind = obj.Interior(XPf);
            
            %fine grid in domain
            XPf = XPf(grid_sq_ind,:);
            
            scatter(XPf(:,1),XPf(:,2),'red');
        end
    end
    
end

