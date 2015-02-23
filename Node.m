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
        function thisNode = Node(numIn, locIn, boundaryIn, ...
                iVal, iValD, iValDD)
            
            thisNode.num = numIn;
            
            thisNode.loc = locIn;
            
            thisNode.onBoundary = boundaryIn;
            
            thisNode.yNew = iVal;
            thisNode.ydNew = iValD;
            thisNode.yddNew = iValDD;
            
            thisNode.yOld = iVal;
            thisNode.ydOld = iValD;
            thisNode.yddOld = iValDD;
            
        end
    end
    
end

