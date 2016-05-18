classdef SystemEQ < handle
    % A class to build and organize systems of equations using FEM
    
    properties
        nNodes
        nbc
        nEle
        nodes
        ele
        bcs
        bigM
        bigC
        bigK
        force
        dynamic_force
        force_const
    end
    
    methods
        %Constructor
        function obj = SystemEQ(nodeSet,elements)
            obj.nNodes = length(nodeSet);
            obj.nbc = 0;
            obj.nEle = length(elements);
            obj.nodes = nodeSet;
            obj.ele = elements(1).empty;
            obj.bcs = BoundaryCondition.empty;
            obj.bigM = zeros(obj.nNodes,obj.nNodes);
            obj.bigC = zeros(obj.nNodes,obj.nNodes);
            obj.bigK = zeros(obj.nNodes,obj.nNodes);
            obj.force = zeros(obj.nNodes,1);
            obj.dynamic_force = false;
            %assemble element matrices into global system
            for i = 1:obj.nEle
                obj.addElement(elements(i));
            end
            
        end
        
        %Add local stiffness matrices to global system
        function [] = addElement(obj,eleToAdd)
            obj.ele(eleToAdd.num) = eleToAdd;
            
            looplim = eleToAdd.nnpe;
            
            %get node numbers for element to be assembled
            nn = zeros(eleToAdd.nnpe,1);
            for i = 1:looplim
                nn(i) = eleToAdd.nodes(i).num;
            end
            
            addK = eleToAdd.eleK;
            addC = eleToAdd.eleC;
            addM = eleToAdd.eleM;
            for k = 1:looplim
                for j = 1:looplim
                    obj.bigK(nn(k),nn(j)) = ...
                        obj.bigK(nn(k),nn(j)) + addK(k,j);
                    obj.bigC(nn(k),nn(j)) = ...
                        obj.bigC(nn(k),nn(j)) + addC(k,j);
                    obj.bigM(nn(k),nn(j)) = ...
                        obj.bigM(nn(k),nn(j)) + addM(k,j);
                end
            end

            eleToAdd.computeForce(0)
            
            %external force assembly
            for m = 1:length(nn)
                obj.force(nn(m)) = obj.force(nn(m)) - eleToAdd.force(m);
            end

            
        end
        
        function [] = updateForce(obj,t)
            
            if obj.dynamic_force == true
                
                obj.force = zeros(length(obj.force),1);
                
                for i = 1:obj.nEle

                    obj.ele(i).computeForce(t)

                    %get node numbers for element to be assembled
                    nn = zeros(obj.ele(i).nnpe,1);
                    for j = 1:length(nn)
                        nn(j) = obj.ele(i).nodes(j).num;
                    end

                    %external force assembly
                    for m = 1:length(nn)
                        obj.force(nn(m)) = obj.force(nn(m)) - obj.ele(i).force(m);
                    end

                end
            else 
                obj.force = obj.force_const;
            end
            
            
        end
        
        function [] = ready(obj)
            
            if obj.dynamic_force == false
                obj.force_const = obj.force;
            end
            
        end
        
        function [] = addBC(obj,BC)
            
            obj.nbc = obj.nbc + 1;
            obj.bcs(BC.num) = BC;
            
        end
            
        function [] = computeDirs(obj)
            
            for i = 1:obj.nEle
                
                obj.ele(i).computeDir()
                
            end
        
        end
        
        function [] = reset(obj)
            for i = 1:obj.nNodes
                obj.nodes(i).y = obj.nodes(i).yic;
                obj.nodes(i).yd = obj.nodes(i).ydic;
                obj.nodes(i).ydd = 0;
            end
            for j = 1:obj.nEle
                obj.ele(j).ypOld = 0;
                obj.ele(j).qOld = 0;
                obj.ele(j).qdOld = 0;
                obj.ele(j).qNew = 0;
                obj.ele(j).qdNew = 0;
            end
            obj.updateForce(0);
        end
        
    end
    
end

