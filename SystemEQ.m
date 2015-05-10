classdef SystemEQ < handle
    %defines class which organizes system of equations for WOOFE
    
    properties
        nNodes = -29;
        nbc = -29;
        nEle = -29;
        J = zeros(-29,-29);
        RHS = zeros(-29,1);
        
        gs4 = GS4(-29,-29,-29,-29,-29);
        
        nodes = cell(2,1);
        ele = cell(2,1);
        bcArray = cell(1,1);
    end
    
    methods
        %Constructor
        function this = SystemEQ(nNo,nEleIn,nodeSet,gs4In)
           
            this.J = zeros(nNo,nNo);
            this.RHS = zeros(nNo,1);
            this.nNodes = nNo;
            this.gs4 = gs4In;
            this.nbc = 0;
            this.nEle = nEleIn;
            
            for i = 1:length(nodeSet)
                this.nodes{i} = nodeSet(i);
            end
            
        end
        
        %Add local stiffness matrices to global system
        function [] = addElement(theSystem,eleToAdd)
                        
            theSystem.ele{eleToAdd.num} = eleToAdd;
            
            looplim = eleToAdd.nNodePerEle;
            
            numDynAdd = eleToAdd.numDyn;
            
            %get node numbers for element to be assembled
            nn = zeros(eleToAdd.nNodePerEle,1);
            for i = 1:looplim
                nn(i) = eleToAdd.nodes{i}.num;
            end
            
            for n = 1:numDynAdd % loop over each elemental addition to J/RHS
                
                %find nodal contributions to J
                addJ = eleToAdd.toJ{n};
                
                for k = 1:looplim
                    for j = 1:looplim
                        theSystem.J(nn(k),nn(j)) = ...
                            theSystem.J(nn(k),nn(j)) + addJ(k,j);
                    end
                end
                
                %find nodal contribution to RHS
                addRHS = eleToAdd.toRHS{n};
                
                for m = 1:looplim
                    theSystem.RHS(nn(m)) = theSystem.RHS(nn(m)) + addRHS(m);
                end
                
            end
                        
            upF = eleToAdd.toRHS{n+1};
            %external force assembly
            for m = 1:looplim
                theSystem.RHS(nn(m)) = theSystem.RHS(nn(m)) - upF(m);
            end
            
        end
        
        function [] = addBC(sys,BC)
            
            sys.nbc = sys.nbc + 1;
            
            sys.bcArray{BC.num} = BC;
            
        end
        
        function [] = changeBC(sys,BC)
            
            sys.bcArray{BC.num} = BC;
            sys.setBCs();
            
        end
        
        
        function [] = setBCs(sys)
            
            for i = 1:sys.nbc
                
                bcToSet = sys.bcArray{i};
                nNo = bcToSet.where;
                rhsVal = bcToSet.changeRHS;
                kVal = bcToSet.addK;
                cVal = bcToSet.addC;
                
                if bcToSet.type == 1
                    
                    sys.J(nNo,:) = zeros(1,length(sys.J(nNo,:)));
                    sys.J(:,nNo) = zeros(length(sys.J(nNo,:)),1);
                    sys.J(nNo,nNo) = 1;
                                        
                    sys.RHS(nNo) = sys.nodes{nNo}.yNew - rhsVal;
                    
                elseif bcToSet.type == 2
                    
                    sys.RHS(nNo) = sys.RHS(nNo) - rhsVal;
                    
                elseif bcToSet.type == 3
                    
                    sys.J(nNo,nNo) = sys.J(nNo,nNo) + kVal;
                    sys.RHS(nNo) = sys.RHS(nNo) + kVal*sys.nodes{nNo}.yNew;
                    sys.RHS(nNo) = sys.RHS(nNo) - rhsVal;
                    
                elseif bcToSet.type == 4
                    
                    sys.J(nNo,nNo) = sys.J(nNo,nNo) + kVal;
                    sys.RHS(nNo) = sys.RHS(nNo) + kVal*sys.nodes{nNo}.yNew;
                    sys.RHS(nNo) = sys.RHS(nNo) + cVal*sys.nodes{nNo}.ydNew;
                    sys.RHS(nNo) = sys.RHS(nNo) - rhsVal;
                    
                end
                
            end
            
        end
        
        function [] = updateRHS(sys)
            
            sys.RHS = zeros(size(sys.RHS));
            
            % update RHS
            for q = 1:sys.nEle
                sys.ele{q}.updateEleRHS(sys.gs4);
                                
                looplim = sys.ele{q}.nNodePerEle;
                
                numDynAdd = sys.ele{q}.numDyn;
                
                %get node numbers for element to be assembled
                nn = zeros(sys.ele{q}.nNodePerEle,1);
                for i = 1:looplim
                    nn(i) = sys.ele{q}.nodes{i}.num;
                end
                
                for n = 1:numDynAdd % loop over each elemental addition to J/RHS
                    
                    %find nodal contribution to RHS
                    addRHS = sys.ele{q}.toRHS{n};
                    
                    for m = 1:looplim
                        sys.RHS(nn(m)) = sys.RHS(nn(m)) + addRHS(m);
                    end
                    
                end
                
                upF = sys.ele{q}.toRHS{n+1};
                %external force assembly
                for m = 1:looplim
                    sys.RHS(nn(m)) = sys.RHS(nn(m)) - upF(m);
                end
                
            end
            
            sys.setBCs()
            
        end
        
        
        function [] = updateSystem(sys)
            
            sys.RHS = zeros(size(sys.RHS));
            sys.J = zeros(size(sys.J));            
            
            for j = 1:sys.nEle
                sys.ele{j}.updateAll(sys.gs4) 
                sys.addElement(sys.ele{j})
            end
            
            sys.setBCs()
            
        end
        
        
        function [] = solve(sys, tol, maxIter)
                        
            delta = 29;
            nlcount = 0;
            
            while norm(delta) > tol
                                         
                delta = -sys.J\sys.RHS;
                                              
                for i = 1:length(delta)
                    sys.nodes{i}.yNew = sys.nodes{i}.yNew + delta(i);
                end
                
                sys.updateRHS()
                %sys.updateSystem()
                
                nlcount = nlcount+1;
                
                if nlcount > maxIter
                    fprintf('\n\n EXCEEDED NL ITER LIMIT\n\n')
                    fprintf('norm(delta) = %e\n',norm(delta))
                    break
                end
                
            end
            
            
            
        end

        function [] = timeMarch(sys)
            
            for i = 1:sys.nNodes
               [newDD,newD] = sys.gs4.march(sys.nodes{i});
               sys.nodes{i}.yddNew = newDD;
               sys.nodes{i}.ydNew = newD;
               
            end
            
            for i = 1:sys.nNodes
                sys.nodes{i}.yOld = sys.nodes{i}.yNew;
                sys.nodes{i}.ydOld = sys.nodes{i}.ydNew;
                sys.nodes{i}.yddOld = sys.nodes{i}.yddNew;
            end
            
            sys.updateSystem()
            
        end
        
        function [] = computeDirs(sys)
            
            for i = 1:sys.nEle
                sys.ele{i}.computeDir()
            end
        
        end    
            
            
    end
    
end

