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
        n
        tn
        tnp1
        tnpw1
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
            
            this.n = 0;
            this.tn = 0;
            this.tnp1 = this.dt;
        end

        
        function [] = time_march(gs4, sys)
          
            delt = gs4.dt;
            
            gs4.tn = delt*gs4.n;
            gs4.n = gs4.n + 1;
            gs4.tnp1 = delt*(gs4.n);
            
            yddn = zeros(sys.nNodes,1);
            ydn = zeros(sys.nNodes,1);
            yn = zeros(sys.nNodes,1);
            yddnp1 = zeros(sys.nNodes,1);
            ydnp1 = zeros(sys.nNodes,1);
            ynp1 = zeros(sys.nNodes,1);
            for i = 1:sys.nNodes
                yddn(i) = sys.nodes(i).yddOld;
                ydn(i) = sys.nodes(i).ydOld;
                yn(i) = sys.nodes(i).yOld;
            end
            
            %make some space for temporary stiffness matrices
            Emm = sys.bigM;
            Cee = sys.bigC;
            Kay = sys.bigK;
            
            % set type 3 BC values before shifting
            for i = 1:sys.nbc
                bc = sys.bcs(i);
                no = bc.where;
                kVal = bc.addK;
                cVal = bc.addC;
                if bc.type == 3
                    Kay(no,no) = Kay(no,no) + kVal;
                    Cee(no,no) = Cee(no,no) + cVal;
                end
            end
            
            
            %perform shift if first order
            if gs4.order == 1
                Emm = Cee;
                Cee = Kay;
                Kay = zeros(sys.nNodes,sys.nNodes);
                yddn = ydn;
                ydn = yn;
            end
            
            LHS = gs4.lam6w1*Emm ...
                + gs4.lam5w2*gs4.dt*Cee + gs4.lam3w3*gs4.dt^2*Kay;
            
            %TODO include time dependent force terms!
            %check if force term is time-depenedent then update as
            %necessary
            
            sys.updateForce();
            
            RHS = - Emm*yddn ...
                  - Cee*(ydn + gs4.lam4w1*delt*yddn) ...
                  - Kay*(yn + gs4.lam1w1*delt*ydn ...
                                + gs4.lam2w2*delt^2*yddn) ...
                                - sys.force;
                  %+ (1-gs4.w1)*sys.fOld + gs4.w1*sys.fNew;
            
            % then set the dirichlet conditions
            for i = 1:sys.nbc
                bc = sys.bcs(i);
                no = bc.where;
                rhsVal = bc.changeRHS;
                RHS(no) = RHS(no) + rhsVal;
                if bc.type == 1
                    LHS(no,:) = zeros(1,length(LHS(no,:)));
                    LHS(:,no) = zeros(length(LHS(no,:)),1);
                    LHS(no,no) = 1;      
                    if gs4.order == 2
                        RHS(no) = ( bc.changeRHS - yn(no) ...
                            - gs4.lam1*gs4.dt*ydn(no) ...
                            - gs4.lam2*gs4.dt^2*yddn(no) ) ...
                            / gs4.lam3/gs4.dt^2;
                    elseif gs4.order == 1
                        RHS(no) = (bc.changeRHS - ydn(no) - delt*yddn(no))...
                                /gs4.lam5/delt;
                    end
                end
            end
            
            %solve!
            delydd = LHS \ RHS;
            
            %compute updates
            for i = 1:sys.nNodes
                yddnp1(i) = yddn(i) + delydd(i);
                ydnp1(i) = ydn(i) + gs4.lam4*delt*yddn(i) ...
                        + gs4.lam5*delt*delydd(i);
                ynp1(i) = yn(i) + gs4.lam1*delt*ydn(i) ...
                        + gs4.lam2*delt^2*yddn(i) ...
                        + gs4.lam3*delt^2*delydd(i);
            end

            %shift back before updating global system
            if gs4.order == 1
                ynp1 = ydnp1;
                ydnp1 = yddnp1;
            end

            for i = 1:sys.nNodes
                sys.nodes(i).yddNew = yddnp1(i);
                sys.nodes(i).ydNew = ydnp1(i);
                sys.nodes(i).yNew = ynp1(i);
                sys.nodes(i).yddOld = sys.nodes(i).yddNew;
                sys.nodes(i).ydOld = sys.nodes(i).ydNew;
                sys.nodes(i).yOld = sys.nodes(i).yNew;
            end
            
            
        end
        
    end
    
end

