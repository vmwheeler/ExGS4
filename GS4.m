classdef GS4 < handle
    %Generalized Single Step Single Solve class definition
    
    properties
        rhoMax = -29;
        rhoMin = -29;
        rhoSpu = -29;
        shift = -29;
        order = -29;
    end
    
    properties (SetAccess = private)
        % algorithm parameters
        lam1w1 = -29;
        lam2w2 = -29;
        lam3w3 = -29;
        lam4w1 = -29;
        lam5w2 = -29;
        lam6w1 = -29;
        w1 = -29;
        lam1 = -29;
        lam2 = -29;
        lam3 = -29;
        lam4 = -29;
        lam5 = -29;
        
        dt = -29;
    end
    
    methods
        %Constructor
        function this = GS4(rmx,rmn,rsp,dtIn,sysOrder)
            
            this.order = sysOrder;
            
            this.dt = dtIn;
            
            if sysOrder == 1
                rmx = 1;
            end
            this.rhoMax = rmx;
            this.rhoMin = rmn;
            this.rhoSpu = rsp;
            
            this.shift =...
                (2+rmn+rmx+rsp-rmn*rmx*rsp)/...
                ((1+rmn)*(1+rmx)*(1+rsp))...
                - (3+rmn+rmx-rmn*rmx)/(2*(1+rmn)*(1+rmx));
                                
            this.lam1w1 = (3+rmn+rmx-rmn*rmx)/(2*(1+rmn)*(1+rmx));
            this.lam2w2 = 1/((1+rmn)*(1+rmx));
            this.lam3w3 = 1/((1+rmn)*(1+rmx)*(1+rsp));
            this.lam4w1 = (3+rmn+rmx-rmn*rmx)/(2*(1+rmn)*(1+rmx));
            this.lam5w2 = 2/((1+rmn)*(1+rmx)*(1+rsp));
            this.lam6w1 = (2+rmn+rmx+rsp-rmn*rmx*rsp)/((1+rmn)*(1+rmx)*(1+rsp));
            this.w1 = (3+rmn+rmx-rmn*rmx)/(2*(1+rmn)*(1+rmx));
            this.lam1 = 1;
            this.lam2 = 1/2;
            this.lam3 = 1/(2*(1+rsp));
            this.lam4 = 1;
            this.lam5 = 1/(1+rsp);
            
        end
        
        function [newDD, newD] = march(gs4,nodeIn)
            
            if gs4.order == 1
                newDD = 0;
                newD = (1-1/gs4.lam5)*nodeIn.ydOld ...
                    + 1/gs4.lam5/gs4.dt*(nodeIn.yNew-nodeIn.yOld);
                
            elseif gs4.order == 2
                
                newDD = (1-gs4.lam2/gs4.lam3)*nodeIn.yddOld ...
                    - 1/gs4.lam3/gs4.dt*nodeIn.ydOld ...
                    + 1/gs4.lam3/gs4.dt/gs4.dt*(nodeIn.yNew-nodeIn.yOld);
                
                newD = (1-gs4.lam5/gs4.lam3)*nodeIn.ydOld ...
                    + (1-gs4.lam2*gs4.lam5/gs4.lam3)*gs4.dt*nodeIn.yddOld...
                    + gs4.lam5/gs4.lam3/gs4.dt*(nodeIn.yNew-nodeIn.yOld);
                
            end
            
            
        end
            
        
        
    end
    
end

