classdef Element1DL < handle
    %The 1D linear implicit GS4 element class for WOOFE
    
    properties
        num = -29;
        nodes = cell(2,1);
        const = [-29,-29];
        nnpe = 2;
        dx = -29;
        toLHS = cell(2,1);
        toRHS = cell(2,1);
        force = [-29,-29];
        numDyn = 2;
        yp = -29;
    end
    
    methods
        %
        % Constructor... send in nodes and material properties
        function this = Element1DL(numIn, nodeInfo, constIn, fIn, gs4)
            
            this.num = numIn;
            
            this.nodes{1} = nodeInfo(1);
            this.nodes{2} = nodeInfo(2);
            
            this.force(1) = fIn(1);
            this.force(2) = fIn(2);
            
            this.const = constIn;
            
            this.dx = abs(this.nodes{1}.loc ...
                - this.nodes{2}.loc );

            this.updateAll(gs4)
            
        end
        %
        
        function [] = updateAll(ele, gs4)
            
            ele.updateEleRHS(gs4)
            ele.updateEleLHS(gs4)
            
        end
        
        function [] = updateEleRHS(ele, gs4)
            
            yw1 = ele.nodes{1}.yOld ...
                + gs4.lam5w2/gs4.lam5 ...
                *(ele.nodes{1}.yNew - ele.nodes{1}.yOld);
            yw2 = ele.nodes{2}.yOld ...
                + gs4.lam5w2/gs4.lam5 ...
                *(ele.nodes{2}.yNew - ele.nodes{2}.yOld);
            
            coeff1 = (1-gs4.lam6w1/gs4.lam5);
            coeff2 = gs4.lam6w1/gs4.lam5/gs4.dt;
            
            ydw1 = coeff1*ele.nodes{1}.ydOld ...
                + coeff2*(ele.nodes{1}.yNew - ele.nodes{1}.yOld);
            ydw2 = coeff1*ele.nodes{2}.ydOld ...
                + coeff2*(ele.nodes{2}.yNew - ele.nodes{2}.yOld);
                        
            eleK = ele.const(2)*1/ele.dx * [ 1 -1 ; ...
                                           -1  1 ];  
                                       
            eleC = ele.const(1)*ele.dx/6 * [ 2 1 ; ...
                                            1 2 ];
                                        
            f1 = ele.dx/2*ele.force(1);
            f2 = ele.dx/2*ele.force(2);
            
            ele.toRHS{1} = eleC*[ydw1;ydw2];
            ele.toRHS{2} = eleK*[yw1;yw2];
            ele.toRHS{3} = [f1;f2];
        end
        
        function [] = updateEleLHS(ele,gs4)         
            
            eleK = ele.const(2)*1/ele.dx * [ 1 -1 ; ...
                                           -1  1 ];  
                                       
            eleC = ele.const(1)*ele.dx/6 * [ 2 1 ; ...
                                            1 2 ];
                                        
            ele.toLHS{1} = gs4.lam6w1/gs4.lam5/gs4.dt*eleC;                                     
            ele.toLHS{2} = gs4.lam5w2/gs4.lam5*eleK;
            
        end
        
        function [] = computeDir(ele)
            % define spatial derivative within element
            ele.yp = (ele.nodes{2}.yNew - ele.nodes{1}.yNew)/ele.dx;
            
        end
        
    end
    
end

