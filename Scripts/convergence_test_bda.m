clear;
close all;
clc;
addpath('../Base','../Elements','../ForceTerms/','../Extras');

%% Physical and numerical constants
numEle = 3;
numNodes = numEle+1;
Kn = 1.0;
FT = 0.0;

%% Initialize vectors
problemIC = zeros(numNodes,1);
sol = zeros(numNodes,1);
flux = zeros(numNodes-1,1);


%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

% initialize (set ICs) and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end

%% Create elements from mesh
% set the force term 
fhandle = @(x,t)ballisticTerm(x,t,Kn);

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = CF_1DL(i,eleNodes, [1,1,Kn^2/3,FT], fhandle);
end

%% Set up system of equations
% initiallize
sysEQ = SystemEQ(nodes,ele);
sysEQ.dynamic_force = true;

%% set BCs
BC1 = BoundaryCondition(1,3,1,0,Kn*0.5,0.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0,Kn*0.5,0.0);
sysEQ.addBC(BC2);


%% Solve!
sysEQ.ready()

%% time for the convergence study
rhoMax = 0.29; rhoMin = 0.28; rhoEss = 0.11;
tEnd = 0.15;
nSteps_exact = 1000;
nSteps_set = [64,32,16,8];
whichnode = 2;

gs4 = GS4(rhoMax,rhoMin,rhoEss,1.,1);
gs4.convergence_test(sysEQ,tEnd,nSteps_exact,nSteps_set,whichnode)



