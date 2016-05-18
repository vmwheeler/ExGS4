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
        dynamic_force
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
            this.tnpw1 = this.tn + this.w1*this.tnp1;
            
            this.dynamic_force = false;
        end

        
        function [] = time_march(gs4, sys)
          
            %some matlab-fu to make sure we are integrating a system
            %of the expected order
            if any(any(sys.bigM))
                gs4.order = 2;
            else
                gs4.order = 1;
            end
            
            delt = gs4.dt;
            
            gs4.tn = delt*gs4.n;
            gs4.n = gs4.n + 1;
            gs4.tnp1 = delt*(gs4.n);
            gs4.tnpw1 = gs4.tn + gs4.w1*(gs4.tnp1-gs4.tn);
            
            
            yddn = zeros(sys.nNodes,1);
            ydn = zeros(sys.nNodes,1);
            yn = zeros(sys.nNodes,1);
            for i = 1:sys.nNodes
                yddn(i) = sys.nodes(i).ydd;
                ydn(i) = sys.nodes(i).yd;
                yn(i) = sys.nodes(i).y;
            end
            
            %make some space for temporary stiffness matrices
            Emm = sys.bigM;
            Cee = sys.bigC;
            Kay = sys.bigK;
            
            sys.updateForce(gs4.tnpw1);
            
            % set type 3 BC values in stiffness matrices before shifting
            for i = 1:sys.nbc
                bc = sys.bcs(i);
                no = bc.where;
                if bc.type == 3
                    sys.force(no) = sys.force(no) - bc.changeRHS;
                    Kay(no,no) = Kay(no,no) + bc.addK;
                    Cee(no,no) = Cee(no,no) + bc.addC;
                elseif bc.type == 2
                    sys.force(no) = sys.force(no) - bc.changeRHS;
                end
            end
            
            %perform shift if first order
            if gs4.order == 1
                Emm = Cee;
                Cee = Kay;
                Kay = zeros(sys.nNodes,sys.nNodes);
                yddn = ydn;
                ydn = yn;
                yn = zeros(sys.nNodes,1);
            end
            

%             
            % enforce a consistent IC for the highest time derivative
            if gs4.tn == 0
                yddn = Emm \ (-sys.force-Cee*ydn-Kay*yn);
            end

            
            LHS = gs4.lam6w1*Emm + gs4.lam5w2*delt*Cee ...
                    + gs4.lam3w3*delt^2*Kay;
            
                
            RHS = - Emm*yddn ...
                  - Cee*(ydn + gs4.lam4w1*delt*yddn) ...
                  - Kay*(yn + gs4.lam1w1*delt*ydn ...
                                + gs4.lam2w2*delt^2*yddn) ...
                                - sys.force;
                       
                        
            
                            
            % then set the dirichlet conditions and add rhs bc terms
            for i = 1:sys.nbc
                bc = sys.bcs(i);
                no = bc.where;
                if bc.type == 1
                    LHS(no,:) = zeros(1,length(LHS(no,:)));
                    LHS(no,no) = 1;      
                    if gs4.order == 2
                        RHS(no) = ( bc.changeRHS - yn(no) ...
                            - gs4.lam1*delt*ydn(no) ...
                            - gs4.lam2*delt^2*yddn(no) ) ...
                            / gs4.lam3/delt^2;
                    elseif gs4.order == 1
                        RHS(no) = (bc.changeRHS - ydn(no) - delt*yddn(no))...
                                /gs4.lam5/delt;
                    end
                end
            end
            
   

            %solve!
            delydd = LHS \ RHS;
            
            %compute updates
            yddnp1 = yddn + delydd;
            ydnp1 = ydn + gs4.lam4*delt*yddn + gs4.lam5*delt*delydd;
            ynp1 = yn + gs4.lam1*delt*ydn + gs4.lam2*delt^2*yddn ...
                    + gs4.lam3*delt^2*delydd;

            %shift back before updating global system
            if gs4.order == 1
                ynp1 = ydnp1;
                ydnp1 = yddnp1;
            end

            for i = 1:sys.nNodes
                sys.nodes(i).ydd = yddnp1(i);
                sys.nodes(i).yd = ydnp1(i);
                sys.nodes(i).y = ynp1(i);
            end
            
            
        end
        
        function [] = convergence_test(obj,sys,tEnd,nEx,nSet,whichnode)
            
            nodes = sys.nodes;
            
            
            dt_exact = tEnd/nEx;
            obj.dt = dt_exact;

            fprintf('Calculating exact solution for case without time shift\n')
            for i = 1:nEx
                obj.time_march(sys);
            end

            tEnpw1 = obj.tnpw1;
            phi = obj.shift;

            yExact = nodes(whichnode).y;
            ydExact = nodes(whichnode).yd;
            yddExact = nodes(whichnode).ydd;

            
            %reset sol to initial state
            sys.reset()

            fprintf('Performing convergence study for case without time shift\n')
            dts = tEnd./nSet;
            err = zeros(length(nSet),1);
            errd = zeros(length(nSet),1);
            errdd = zeros(length(nSet),1);
            for j = 1:length(nSet)
                %first non-shifted
                obj.dt = dts(j);
                obj.reset();
                sys.ready();
                for i = 1:nSet(j)
                    obj.time_march(sys);
                end
                
                err(j) = abs(nodes(whichnode).y - yExact);
                errd(j) = abs(nodes(whichnode).yd - ydExact);
                errdd(j) = abs(nodes(whichnode).ydd - yddExact);
                sys.reset()
            end

                        
            %now for shifted study

            fprintf('Calculating exact solution for case with time shift\n')
            obj.dt = tEnd/(nEx-phi);
            sys.reset()
            obj.reset()
            sys.ready()
            for i = 1:nEx
                obj.time_march(sys);
            end

            yExact_sh = nodes(whichnode).y;
            ydExact_sh = nodes(whichnode).yd;
            yddExact_sh = nodes(whichnode).ydd;

            sys.reset()
            fprintf('Performing convergence study for case without time shift\n\n')
            dts_sh = tEnd./(nSet-phi);
            err_sh = zeros(length(nSet),1);
            errd_sh = zeros(length(nSet),1);
            errdd_sh = zeros(length(nSet),1);
            for j = 1:length(nSet)
                %then shifted
                obj.dt = dts_sh(j);
                obj.reset();
                %gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,2);
                sys.ready();
                for i = 1:nSet(j)
                    obj.time_march(sys);
                end
                err_sh(j) = abs(nodes(whichnode).y - yExact_sh);
                errd_sh(j) = abs(nodes(whichnode).yd - ydExact_sh);
                errdd_sh(j) = abs(nodes(whichnode).ydd - yddExact_sh);
                sys.reset()
            end
            
            %fit line to dt vs error data to get convergence rate (slope)
            ylinfit = polyfit(log(dts),log(err'),1);
            yrate = ylinfit(1);
            yfitvals = polyval(ylinfit,log(dts));

            ydlinfit = polyfit(log(dts),log(errd'),1);
            ydrate = ydlinfit(1);
            ydfitvals = polyval(ydlinfit,log(dts));

            yddlinfit = polyfit(log(dts),log(errdd'),1);
            yddrate = yddlinfit(1);
            yddfitvals = polyval(yddlinfit,log(dts));
            
            %same for shifted values
            ylinfit_sh = polyfit(log(dts_sh),log(err_sh'),1);
            yrate_sh = ylinfit_sh(1);
            yfitvals_sh = polyval(ylinfit_sh,log(dts_sh));

            ydlinfit_sh = polyfit(log(dts_sh),log(errd_sh'),1);
            ydrate_sh = ydlinfit_sh(1);
            ydfitvals_sh = polyval(ydlinfit_sh,log(dts_sh));

            yddlinfit_sh = polyfit(log(dts_sh),log(errdd_sh'),1);
            yddrate_sh = yddlinfit_sh(1);
            yddfitvals_sh = polyval(yddlinfit_sh,log(dts_sh));

            fprintf('*******rates with no time shift*******\n')
            fprintf('theta: %f\n',yrate)
            fprintf('dtheta: %f\n',ydrate)
            fprintf('ddtheta: %f\n\n',yddrate)
            fprintf('*******rates with time shift*******\n')
            fprintf('theta: %f\n',yrate_sh)
            fprintf('dtheta: %f\n',ydrate_sh)
            fprintf('ddtheta: %f\n\n',yddrate_sh)
            
            %make some figures and display the rates in the legend
            figure('Name','Convergence plot with no time shift')
            set(gca,'fontsize',8)
            set(gcf, 'Units','centimeters', 'Position',[0 0 6 5])
            cy = loglog(dts',err,'vr');
            hold on
            cyd = loglog(dts',errd,'sr');
            cydd = loglog(dts',errdd,'or');
            fy = loglog(dts',exp(yfitvals),'k');
            fyd = loglog(dts',exp(ydfitvals),'k');
            fydd = loglog(dts',exp(yddfitvals),'k');
%             le0 = strcat('$\theta:\;$ ', num2str(yrate));
%             le1 = strcat('$\dot{\theta}:\;$ ', num2str(ydrate));
%             le2 = strcat('$\ddot{\theta}:\;$ ', num2str(yddrate));
            le0 = '$\theta$';
            le1 = '$\dot{\theta}$';
            le2 = '$\ddot{\theta}$';
            legend([cy,cyd,cydd,fy],...
                {le0,le1,le2,'fit'},...
                'interpreter','latex','fontsize',12,'Location','Best')
            
            
            figure('Name','Convergence plot with time shift')
            set(gca,'fontsize',8)
            set(gcf, 'Units','centimeters', 'Position',[0 0 6 5])
            cy_sh = loglog(dts_sh',err_sh,'vr');
            hold on
            cyd_sh = loglog(dts_sh',errd_sh,'sr');
            cydd_sh = loglog(dts_sh',errdd_sh,'or');
            fy_sh = loglog(dts_sh',exp(yfitvals_sh),'k');
            fyd_sh = loglog(dts_sh',exp(ydfitvals_sh),'k');
            fydd_sh = loglog(dts_sh',exp(yddfitvals_sh),'k');
%             le0 = strcat('$\theta:\;$ ', num2str(yrate_sh));
%             le1 = strcat('$\dot{\theta}:\;$ ', num2str(ydrate_sh));
%             le2 = strcat('$\ddot{\theta}:\;$ ', num2str(yddrate_sh));
            le0 = '$\theta$';
            le1 = '$\dot{\theta}$';
            le2 = '$\ddot{\theta}$';
            legend([cy_sh,cyd_sh,cydd_sh,fy_sh],...
                {le0,le1,le2,'fit'},...
                'interpreter','latex','fontsize',12,'Location','Best')
            
        end
        
        function [] = reset(this)
            this.n = 0;
            this.tn = 0;
            this.tnp1 = this.dt;
            this.tnpw1 = this.tn + this.w1*this.tnp1;
        end
        
    end
    
end

