classdef BoundaryCondition < handle
    % an object that keeps track of boundary conditions and 
    
    properties
        num = -29
        type = -29;
        where = -29;
        addC = -29;
        addK = -29;
        changeRHS = -29
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

