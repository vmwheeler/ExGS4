classdef Element1DQ < handle
    %The 1D quadratic element class for WOOFE
    
    properties
        num = -29;
        nodes = cell(3,1);
        const = [-29,-29];
        nNodePerEle = 3;
        dx = -29;
        toJ = cell(3,1);
        toRHS = cell(3,1);
        force = [-29,-29,-29];
        numDyn = 2;
        yp = -29;
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        %
        % Constructor... send in nodes and material properties
        function thisEle = Element1DQ(numIn, nodeInfo, constIn, fIn, gs4)
            
            thisEle.num = numIn;
            
            thisEle.nodes{1} = nodeInfo(1);
            thisEle.nodes{2} = nodeInfo(2);
            thisEle.nodes{3} = nodeInfo(3);
            
            thisEle.force(1) = fIn(1);
            thisEle.force(2) = fIn(2);
            thisEle.force(3) = fIn(3);
            
            thisEle.const = constIn;
            
            thisEle.dx = abs(thisEle.nodes{1}.loc ...
                - thisEle.nodes{3}.loc );

            thisEle.updateAll(gs4)
            
        end
        %
        
        function [] = updateAll(ele, gs4)
            
            ele.updateEleRHS(gs4)
            ele.updateEleJ(gs4)
            
        end
        
        function [] = updateEleRHS(ele, gs4)
            
            coeff1 = (1-gs4.lam6w1/gs4.lam5);
            
            coeff2 = gs4.lam6w1/gs4.lam5/gs4.dt;
            
            
            for i = 1:ele.nNodePerEle
                yw(i) = ele.nodes{i}.yOld ...
                    + gs4.lam5w2/gs4.lam5 ...
                    *(ele.nodes{i}.yNew - ele.nodes{i}.yOld);
            
                ydw(i) = coeff1*ele.nodes{i}.ydOld ...
                + coeff2*(ele.nodes{i}.yNew - ele.nodes{i}.yOld);
            end
                                       
            eleC = ele.const(1)*ele.dx/30 * [ 4,  2, -1; ...
                                              2,  16, 2;...
                                             -1,  2,  4];
                                        
            eleK = ele.const(2)*1/3/ele.dx * [ 7,  -8,  1; ...
                                              -8,  16, -8;...
                                               1,  -8,  7];                              
                                        
            f1 = ele.dx/6*ele.force(1);
            f2 = 2*ele.dx/3*ele.force(2);
            f3 = ele.dx/6*ele.force(3);
            
            ele.toRHS{1} = eleC*[ydw(1);ydw(2);ydw(3)];
            ele.toRHS{2} = eleK*[yw(1);yw(2);yw(3)];
            ele.toRHS{3} = [f1;f2;f3];
        end
        
        function [] = updateEleJ(ele,gs4)         
            
            eleC = ele.const(1)*ele.dx/30 * [ 4,  2, -1; ...
                                              2,  16, 2;...
                                             -1,  2,  4];
                                        
            eleK = ele.const(2)*1/3/ele.dx * [ 7,  -8,  1; ...
                                              -8,  16, -8;...
                                               1,  -8,  7]; 
                                        
            ele.toJ{1} = gs4.lam6w1/gs4.lam5/gs4.dt*eleC;                                     
            ele.toJ{2} = gs4.lam5w2/gs4.lam5*eleK;
            
        end
        
        function [] = computeDir(ele)
            % define spatial derivative within element
            ele.yp = (ele.nodes{3}.yNew - ele.nodes{1}.yNew)/ele.dx;
            
        end
        
    end
    
end

