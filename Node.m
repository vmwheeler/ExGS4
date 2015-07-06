classdef Node < handle
    %The general node class
    
    properties
        loc
        num
        onBoundary
        yNew
        ydNew
        yddNew
        yOld
        ydOld
        yddOld
    end
    
    methods
        function this = Node(numIn, locIn, boundaryIn, ic, icD, icDD)
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

