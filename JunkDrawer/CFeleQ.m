classdef CFeleQ < handle
    %The 1D quadratic CF element class for WOOFE
    
    properties
        num = -29;
        nodes = cell(3,1);
        const = [-29,-29,-29,-29];
        nNodePerEle = 3;
        dx = -29;
        toJ = cell(4,1);
        toRHS = cell(4,1);
        force = [-29,-29,-29];
        numDyn = 3;
        yp = -29;
        ydp = -29;
    end
    
    properties (SetAccess = private)
        
    end
    
    methods
        %
        % Constructor... send in nodes and material properties
        function thisEle = CFeleQ(numIn, nodeInfo, constIn, fIn, gs4)
            
            thisEle.num = numIn;
            
            thisEle.nodes{1} = nodeInfo(1);
            thisEle.nodes{2} = nodeInfo(2);
            thisEle.nodes{3} = nodeInfo(3);
            
            thisEle.force(1) = fIn(1);
            thisEle.force(2) = fIn(2);
            thisEle.force(3) = fIn(2);
            
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
            
            cc = gs4.lam3w3/gs4.lam3;
            ccd = gs4.lam5w2/gs4.lam3;
            ccdd = gs4.lam6w1/gs4.lam3;

            
            for i = 1:ele.nNodePerEle
                yw(i) = ele.nodes{i}.yOld ...
                    + gs4.dt^2*(gs4.lam2w2 - gs4.lam2*cc)*ele.nodes{i}.yddOld...
                    + (gs4.lam1w1  - cc)*gs4.dt*ele.nodes{i}.ydOld ...
                    + cc*(ele.nodes{i}.yNew - ele.nodes{i}.yOld);
                
                ydw(i) = (1-ccd)*ele.nodes{i}.ydOld ...
                    + (gs4.lam4w1-gs4.lam2*ccd)*gs4.dt*ele.nodes{i}.yddOld...
                    + ccd/gs4.dt*(ele.nodes{i}.yNew - ele.nodes{i}.yOld);
                
                yddw(i) = (1-gs4.lam2*ccdd)*ele.nodes{i}.yddOld ...
                    - ccdd/gs4.dt*ele.nodes{i}.ydOld ...
                    + ccdd/gs4.dt^2*(ele.nodes{i}.yNew - ele.nodes{i}.yOld);
            end
            
            
            eleM = ele.const(1)*ele.dx/30 * [ 4,  2, -1; ...
                                              2,  16, 2;...
                                             -1,  2,  4];  
                                       
            eleC = ele.const(2)*ele.dx/30 * [ 4,  2, -1; ...
                                              2,  16, 2;...
                                             -1,  2,  4];
                                        
            eleK = ele.const(3)*1/3/ele.dx * [ 7,  -8,  1; ...
                                              -8,  16, -8;...
                                               1,  -8,  7]; 
                                           
            eleC2 = ele.const(3)*ele.const(4)*1/3/ele.dx * [ 7,  -8,  1; ...
                                              -8,  16, -8;...
                                               1,  -8,  7];                               
                                        
            f1 = ele.dx/6*ele.force(1);
            f2 = 2*ele.dx/3*ele.force(2);
            f3 = ele.dx/6*ele.force(3);
            
            ele.toRHS{1} = (eleC+eleC2)*[ydw(1);ydw(2),;ydw(3)];
            ele.toRHS{2} = eleK*[yw(1);yw(2),;yw(3)];
            ele.toRHS{3} = eleM*[yddw(1);yddw(2);yddw(3)];
            ele.toRHS{4} = [f1;f2;f3];
        end
        
        function [] = updateEleJ(ele,gs4)         
            
            eleM = ele.const(1)*ele.dx/30 * [ 4,  2, -1; ...
                                              2,  16, 2;...
                                             -1,  2,  4];  
                                       
            eleC = ele.const(2)*ele.dx/30 * [ 4,  2, -1; ...
                                              2,  16, 2;...
                                             -1,  2,  4];
                                        
            eleK = ele.const(3)*1/3/ele.dx * [ 7,  -8,  1; ...
                                              -8,  16, -8;...
                                               1,  -8,  7];    
                                         
            eleC2 = ele.const(3)*ele.const(4)*1/3/ele.dx * [ 7,  -8,  1; ...
                                              -8,  16, -8;...
                                               1,  -8,  7];
                                        
            ele.toJ{1} = gs4.lam6w1/gs4.lam3/gs4.dt^2*eleM;                                     
            ele.toJ{2} = gs4.lam5w2/gs4.lam3/gs4.dt*(eleC+eleC2);
            ele.toJ{3} = gs4.lam3w3/gs4.lam3*eleK;
            
        end
        
        function [] = computeDir(ele)
            % define spatial derivative within element
            
            ele.yp = (ele.nodes{3}.yNew - ele.nodes{1}.yNew)/ele.dx;
            ele.ydp = (ele.nodes{3}.ydNew - ele.nodes{1}.ydNew)/ele.dx;
            
        end
        
        
    end
    
end

