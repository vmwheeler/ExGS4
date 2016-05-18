classdef BoundaryCondition < handle
    % an object that keeps track of boundary conditions and 
    
    properties
        num
        type
        where
        addC
        addK
        changeRHS
    end
    
    methods
        
        function this = BoundaryCondition(...
                    numIn,...
                    typeIn,...
                    nodeLocIn,...
                    addCin,...
                    addKin,...
                    repRHSin)
            
            this.num = numIn;
            this.type = typeIn;
            this.where = nodeLocIn;
            this.addC = addCin;
            this.addK = addKin;
            this.changeRHS = repRHSin;
        
        end
            
    end
    
end

