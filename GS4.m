classdef GS4 < handle
    %Generalized Single Step Single Solve class definition
    
    properties

        lam1w1
        lam2w2
        lam3w3
        lam4w1
        lam5w2
        lam6w1
        w1
        lam1
        lam2
        lam3
        lam4
        lam5
        dt
        order
        shift
        
    end
    
    methods
        %Constructor
        function this = GS4(rmx,rmn,rsp,dtIn,sysOrder)
            
            this.order = sysOrder;
            
            this.dt = dtIn;
            
            if sysOrder == 1
                rmx = 1;
            end
            
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

        
        function [] = time_march(gs4, sys)
          
            delt = gs4.dt;
            
            

            
            %perform shift
            Emm = sys.bigC;
            Cee = sys.bigK;
            Kay = zeros(sys.nNodes,sys.nNodes);
            
            yddn = zeros(sys.nNodes,1);
            ydn = zeros(sys.nNodes,1);
            yn = zeros(sys.nNodes,1);
            for i = 1:sys.nNodes
                yddn(i) = sys.nodes(i).yddOld;
                ydn(i) = sys.nodes(i).ydOld;
                yn(i) = sys.nodes(i).yOld;
            end
            

            LHS = gs4.lam6w1*Emm ...
                + gs4.lam5w2*gs4.dt*Cee + gs4.lam3w3*gs4.dt^2*Kay;
            
            %TODO include force terms!
            
            RHS = - Emm*yddn ...
                  - Cee*(ydn + gs4.lam4w1*delt*yddn) ...
                  - Kay*(yn + gs4.lam1w1*delt*ydn ...
                                + gs4.lam2w2*delt^2*yddn); %...
                  %+ (1-gs4.w1)*sys.fOld + gs4.w1*sys.fNew;
            

            %TODO need to write some code here to enforce BCs
            %should be "straight forward" since sysEQ was passed in
            
            %TODO first set neumann and robin conditions
            
            % then set the dirichlet conditions
            for i = 1:sys.nbc
                bc = sys.bcs(i);
                if bc.type == 1
                    no = bc.where;
                    LHS(no,:) = zeros(1,length(LHS(no,:)));
                    LHS(:,no) = zeros(length(LHS(no,:)),1);
                    LHS(no,no) = 1;                   
                    RHS(no) = (bc.changeRHS - ydn(no) - delt*yddn(no))/gs4.lam5/delt;
                end
            end

            
            delydd = LHS \ RHS;
            
            for i = 1:sys.nNodes
                sys.nodes(i).yddNew = yddn(i) + delydd(i);
                sys.nodes(i).ydNew = ydn(i) + gs4.lam4*delt*yddn(i) ...
                        + gs4.lam5*delt*delydd(i);
                sys.nodes(i).yNew = yn(i) + gs4.lam1*delt*ydn(i) ...
                        + gs4.lam2*delt^2*yddn(i) ...
                        + gs4.lam3*delt^2*delydd(i);
                sys.nodes(i).yddOld = sys.nodes(i).yddNew;
                sys.nodes(i).ydOld = sys.nodes(i).ydNew;
                sys.nodes(i).yOld = sys.nodes(i).yNew;
            end
            
            
        end
        
    end
    
end

