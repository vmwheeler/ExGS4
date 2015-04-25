classdef Node < handle
    %The general node class for WOOFE
    
    properties
        loc = -29; %location in space
        num = -29; %node number
        onBoundary = -29;
        yNew = -29;
        ydNew = -29;
        yddNew = -29;
        yOld = -29;
        ydOld = -29;
        yddOld = -29
    end
    
    methods
        function this = Node(numIn, locIn, boundaryIn, ...
                ic, icD, icDD)
            
            this.num = numIn;
            
            this.loc = locIn;
            
            this.onBoundary = boundaryIn;
            
            this.yNew = ic;
            this.ydNew = icD;
            this.yddNew = icDD;
            
            this.yOld = ic;
            this.ydOld = icD;
            this.yddOld = icDD;
            
        end
    end
    
end

