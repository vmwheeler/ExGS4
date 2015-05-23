classdef Element1DL < handle
    %Class describing a 1-D linear finite element
    
    properties
        num %element number
        nodes
        const
        nnpe
        dx
        toLHS
        toRHS
        force
        numDyn
        yp
        eleC
        eleK
    end
    
    methods
        % Constructor... send in nodes and material properties
        function obj = Element1DL(numIn, nodeInfo, constIn, fIn, gs4)
            obj.num = numIn;
            obj.nodes = nodeInfo;
            obj.const = constIn;
            obj.nnpe = 2;
            obj.dx = abs(obj.nodes(1).loc - obj.nodes(2).loc );
            obj.toLHS = cell(2,1);
            obj.toRHS = cell(2,1);
            obj.force = fIn;
            obj.numDyn = 2;
            obj.yp = 0;
            obj.eleC = obj.const(1)*obj.dx/6 * [ 2 1 ; ...
                                                 1 2 ];
            obj.eleK = obj.const(2)*1/obj.dx * [ 1 -1 ; ...
                                                -1  1 ]; 
            obj.updateAll(gs4)
        end
        
        function [] = updateAll(obj, gs4)
            obj.updateEleRHS(gs4)
            obj.updateEleLHS(gs4)
        end
        
        function [] = updateEleRHS(obj, gs4)
            
            yw1 = obj.nodes(1).yOld ...
                + gs4.lam5w2/gs4.lam5 ...
                *(obj.nodes(1).yNew - obj.nodes(1).yOld);
            yw2 = obj.nodes(2).yOld ...
                + gs4.lam5w2/gs4.lam5 ...
                *(obj.nodes(2).yNew - obj.nodes(2).yOld);
            
            coeff1 = (1-gs4.lam6w1/gs4.lam5);
            coeff2 = gs4.lam6w1/gs4.lam5/gs4.dt;
            
            ydw1 = coeff1*obj.nodes(1).ydOld ...
                + coeff2*(obj.nodes(1).yNew - obj.nodes(1).yOld);
            ydw2 = coeff1*obj.nodes(2).ydOld ...
                + coeff2*(obj.nodes(2).yNew - obj.nodes(2).yOld);
                                        
            f1 = obj.dx/2*obj.force(1);
            f2 = obj.dx/2*obj.force(2);
            
            obj.toRHS{1} = obj.eleC*[ydw1;ydw2];
            obj.toRHS{2} = obj.eleK*[yw1;yw2];
            obj.toRHS{3} = [f1;f2];
        end
        
        function [] = updateEleLHS(obj,gs4)         
                                        
            obj.toLHS{1} = gs4.lam6w1/gs4.lam5/gs4.dt*obj.eleC;                                     
            obj.toLHS{2} = gs4.lam5w2/gs4.lam5*obj.eleK;
            
        end
        
        function [] = computeDir(obj)
            % define spatial derivative within element
            obj.yp = (obj.nodes(2).yNew - obj.nodes(1).yNew)/obj.dx;
            %obj.yp = (obj.nodes(2).ydNew - obj.nodes(1).ydNew)/obj.dx;
            
        end
        
    end
    
end

