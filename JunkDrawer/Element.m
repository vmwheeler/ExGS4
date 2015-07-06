classdef Element < handle
    
   properties
        num = 0;
        nodes = Node(0,0,0,0,0,0).empty
        const = [0,0];
        nnpe = 2;
        dx = 0;
        eleM = [0,0;
                0,0]
        eleC = [0,0;
                0,0]
        eleK = [0,0;
                0,0]
        toLHS = cell(2,1);
        toRHS = cell(2,1);
        force = [-29,-29];
        numDyn = 2;
        yp = 0;
    end 
    
end