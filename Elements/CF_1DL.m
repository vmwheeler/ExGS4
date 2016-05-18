classdef CF_1DL < handle
    %Class describing a 1-D linear finite element
    
    properties
        % this is all of the stuff that defines the element
        num 
        nodes
        const
        nnpe
        dx
        force
        fhandle
        yp
        eleM
        eleC
        eleK
        ypOld
        qOld
        qdOld
        qNew
        qdNew
    end
    
    methods
        % Constructor... send in nodes and material properties
        function obj = CF_1DL(numIn, nodeInfo, constIn, fhandle)
            obj.num = numIn;
            obj.nodes = nodeInfo;
            obj.const = constIn;
            obj.nnpe = 2;
            obj.dx = abs(obj.nodes(1).loc - obj.nodes(2).loc );
            obj.yp = 0;
            % a check for a first or second order system
            if obj.const(4) < 1
                obj.eleM = obj.const(1)*obj.dx/6 * [ 2 1 ; ...
                                                     1 2 ];
                eleC2 = obj.const(3)*obj.const(4)*1/obj.dx * [ 1 -1 ; ...
                                                              -1  1 ];
            else
                obj.eleM = [ 0 0 ; ...
                             0 0 ];
                eleC2 = [ 0 0 ; ...
                          0 0 ];
            end
                
            eleC1 = obj.const(2)*obj.dx/6 * [ 2 1 ; ...
                                              1 2 ];
            
            obj.eleC = eleC1 + eleC2;
            obj.eleK = obj.const(3)*1/obj.dx * [ 1 -1 ; ...
                                                -1  1 ];
            
                                            
            obj.fhandle = fhandle;
            obj.force = [0;0];
            
            obj.ypOld = 0;
            obj.qOld = 0; obj.qNew = 0;
            obj.qdOld = 0; obj.qdNew = 0;
 
        end
        
        function [] = view(obj)
            % displays the relevant information about the element
            
            fprintf('--------------------------------------------------------\n')
            fprintf('I am a 1-D linear element with the following properties:\n')
            fprintf('--------------------------------------------------------\n')
            fprintf('Element #: %i\n', obj.num)
            fprintf('Nodes: %i, %i\n', obj.nodes(1).num, obj.nodes(2).num)
            fprintf('M^(e): [ %f, %f ]\n', obj.eleK(1,1), obj.eleK(1,2) )
            fprintf('       [ %f, %f ]\n', obj.eleK(2,1), obj.eleK(2,2) )
            fprintf('K^(e): [ %f, %f ]\n', obj.eleK(1,1), obj.eleK(1,2) )
            fprintf('       [ %f, %f ]\n', obj.eleK(2,1), obj.eleK(2,2) )
            fprintf('C^(e): [ %f, %f ]\n', obj.eleC(1,1), obj.eleC(1,2) )
            fprintf('       [ %f, %f ]\n', obj.eleC(2,1), obj.eleC(2,2) )
            fprintf('y'': %f\n', obj.yp)
            fprintf('Force: [ %e ]\n', obj.force(1) )
            fprintf('       [ %e ]\n', obj.force(2) )
            
        end
        
        function [] = computeForce(obj,t)
            % evaluate integrals to get forcing function contributions to
            % each node
            % note the need for the 'ArrayValued' flag in order to handle a
            % call to integral() in the force term definition

            x1 = obj.nodes(1).loc;
            x2 = obj.nodes(2).loc;
            ff1 = @(x) (1-(x2-x)./obj.dx).*obj.fhandle(x,t);
            ff2 = @(x) (x2-x)./obj.dx.*obj.fhandle(x,t);
            obj.force = [ integral(ff1,x1,x2,'ArrayValued', true);
                          integral(ff2,x1,x2,'ArrayValued', true) ];
            
        end
        
        function [] = computeDir(obj)
            % define spatial derivative within element
            obj.yp = (obj.nodes(2).y - obj.nodes(1).y)/obj.dx;
            
        end
        
        function [] = computeFlux(obj,gs4,Kn,FT)
            % this funtion computes the heat flux for models with a time
            % derivative appearing in their definition

            obj.qdNew = ( -Kn/3*(1-FT)*(obj.ypOld+gs4.w1*(obj.yp-obj.ypOld))...
                     +(gs4.lam6w1+gs4.lam5w2*gs4.dt-gs4.w1*gs4.dt-1)*obj.qdOld ...
                     - obj.qOld)/(gs4.lam6w1+gs4.lam5w2*gs4.dt);
                 
            obj.qNew = obj.qOld + gs4.dt*obj.qdOld ...
                + gs4.lam5*gs4.dt*(obj.qdNew-obj.qdOld); 
            
            obj.qNew = obj.qNew - FT*Kn/3*obj.yp;
            
            obj.qdOld = obj.qdNew;
            obj.qOld = obj.qNew;
            obj.ypOld = obj.yp;
            
        end
            
    end
    
end

