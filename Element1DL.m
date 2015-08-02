classdef Element1DL < handle
    %Class describing a 1-D linear finite element
    
    properties
        num %element number
        nodes
        const
        nnpe
        dx
        force
        fhandle
        yp
        eleC
        eleK
    end
    
    methods
        % Constructor... send in nodes and material properties
        function obj = Element1DL(numIn, nodeInfo, constIn, fhandle)
            obj.num = numIn;
            obj.nodes = nodeInfo;
            obj.const = constIn;
            obj.nnpe = 2;
            obj.dx = abs(obj.nodes(1).loc - obj.nodes(2).loc );
            obj.yp = 0;
            obj.eleC = obj.const(1)*obj.dx/6 * [ 2 1 ; ...
                                                 1 2 ];
            obj.eleK = obj.const(2)*1/obj.dx * [ 1 -1 ; ...
                                                -1  1 ];
            obj.fhandle = fhandle;
            obj.force = [0;0];
 
        end
        
        function [] = computeForce(obj)
            
            % evaluate integrals to get forcing function contributions to
            % each node
            %TODO check if this works for real interesting forcing
            %fucntions... this uses local coordinates so i think it will
            %break down for non-constant forces
            x1 = obj.nodes(1).loc;
            x2 = obj.nodes(2).loc;
            ff1 = @(x) (1-(x2-x)./obj.dx).*obj.fhandle(x);
            ff2 = @(x) ((x2-x)./obj.dx).*obj.fhandle(x);
            obj.force = [ integral(ff1,x1,x2);
                          integral(ff2,x1,x2) ];
            
        end
        
        function [] = computeDir(obj)
            % define spatial derivative within element
            obj.yp = (obj.nodes(2).yNew - obj.nodes(1).yNew)/obj.dx;
            
        end
        
    end
    
end

