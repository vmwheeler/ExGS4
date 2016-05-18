classdef Node < handle
    %The general node class
    
    properties
        loc
        num
        y
        yd
        ydd
        yic
        ydic
    end
    
    methods
        function this = Node(numIn, locIn, ic, icD)
            this.num = numIn;
            this.loc = locIn;
            this.y = ic;
            this.yd = icD;
            this.ydd = 0;
            this.yic = ic;
            this.ydic = icD;
        end
    end
    
end

