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
        bigC
        bigK
        bigM
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
        end
        
        %Add local stiffness matrices to global objtem
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
            for k = 1:looplim
                for j = 1:looplim
                    obj.bigK(nn(k),nn(j)) = ...
                        obj.bigK(nn(k),nn(j)) + addK(k,j);
                    obj.bigC(nn(k),nn(j)) = ...
                        obj.bigC(nn(k),nn(j)) + addC(k,j);
                end
            end
            
            %TODO assemble force
            %external force assembly
            %for m = 1:looplim
            %    obj.RHS(nn(m)) = obj.RHS(nn(m)) - upF(m);
            %end
            
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
        
        
        function [] = setBCs(obj)
            
            for i = 1:obj.nbc
                
                bcToSet = obj.bcArray{i};
                nNo = bcToSet.where;
                rhsVal = bcToSet.changeRHS;
                kVal = bcToSet.addK;
                cVal = bcToSet.addC;
                
                if bcToSet.type == 1
                    
                    obj.LHS(nNo,:) = zeros(1,length(obj.LHS(nNo,:)));
                    obj.LHS(:,nNo) = zeros(length(obj.LHS(nNo,:)),1);
                    obj.LHS(nNo,nNo) = 1;    
                    
                    obj.RHS(nNo) = obj.nodes(nNo).yNew - rhsVal;
                    
                elseif bcToSet.type == 2
                    
                    obj.RHS(nNo) = obj.RHS(nNo) - rhsVal;
                    
                elseif bcToSet.type == 3
                    
                    obj.LHS(nNo,nNo) = obj.LHS(nNo,nNo) + kVal;
                    obj.RHS(nNo) = obj.RHS(nNo) + kVal*obj.nodes(nNo).yNew;
                    obj.RHS(nNo) = obj.RHS(nNo) - rhsVal;
                    
                elseif bcToSet.type == 4 
                    
                    obj.LHS(nNo,nNo) = obj.LHS(nNo,nNo) + kVal;
                    obj.RHS(nNo) = obj.RHS(nNo) + kVal*obj.nodes(nNo).yNew;
                    obj.RHS(nNo) = obj.RHS(nNo) + cVal*obj.nodes(nNo).ydNew;
                    obj.RHS(nNo) = obj.RHS(nNo) - rhsVal;
                    
                end
            end
        end
            
        function [] = computeDirs(obj)
            
            for i = 1:obj.nEle
                
                obj.ele(i).computeDir()
                
            end
        
        end    
            
        function [] = lump(obj)
            
            for i = 1:obj.nNodes
                
                newVal = sum(obj.LHS(i,:))/length(obj.LHS(i,:));
                obj.LHS(i,:) = zeros(1,length(obj.LHS(i,:)));
                obj.LHS(i,i) = newVal;
                
            end
            
        end
        
    end
    
end

