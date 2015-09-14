classdef SystemEQ < handle
    % A class to build and organize systems of equations using FEM
    
    properties
        nNodes
        nbc
        nEle
        nodes
        ele
        bcArray
        bcs
        gs4
        bigM
        bigC
        bigK
        force
    end
    
    methods
        %Constructor
        function obj = SystemEQ(nNo,nEleIn,nodeSet,gs4In,element_type)
            obj.nNodes = nNo;
            obj.nbc = 0;
            obj.nEle = nEleIn;
            obj.gs4 = gs4In;
            obj.nodes = Node.empty;
            obj.ele = element_type.empty;
            for i = 1:length(nodeSet)
                obj.nodes(i) = nodeSet(i);
            end
            obj.bcs = BoundaryCondition.empty;
            obj.bigM = zeros(obj.nNodes,obj.nNodes);
            obj.bigC = zeros(obj.nNodes,obj.nNodes);
            obj.bigK = zeros(obj.nNodes,obj.nNodes);
            obj.force = zeros(obj.nNodes,1);
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
            
            
            
        end
        
        function [] = updateForce(obj)

            obj.force = zeros(length(obj.force),1);
            
            for i = 1:obj.nEle
                
                obj.ele(i).computeForce()
                
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
            
        end
        
        function [] = addBC(obj,BC)
            
            obj.nbc = obj.nbc + 1;
            obj.bcArray{BC.num} = BC;      
            obj.bcs(BC.num) = BC;
            
        end
        
        function [] = changeBC(obj,BC)
            
            obj.bcArray{BC.num} = BC;
            obj.setBCs();
            
        end
            
        function [] = computeDirs(obj)
            
            for i = 1:obj.nEle
                
                obj.ele(i).computeDir()
                
            end
        
        end
        
    end
    
end

