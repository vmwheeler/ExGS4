% A script for performing a convergence test on a chosen problem
clear;
addpath('../Base','../Elements','../ForceTerms/');

%% Physical and numerical constants
numEle = 5;
numNodes = numEle+1;
Kn = 1.0;
FT = 1.0;

%% Initialize vectors
problemIC = zeros(numNodes,1);
sol = zeros(numNodes,1);
flux = zeros(numNodes-1,1);


%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

xPlot = linspace(0,1,numNodes);

problemIC(1) = 1;
% initialize and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end

% set the force term (simply 0 here)
fhandle = @force_term_0;

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = CF_1DL(i,eleNodes, [1,1,Kn^2/3,FT], fhandle);
end

% initiallize
sysEQ = SystemEQ(nodes,ele);

BC1 = BoundaryCondition(1,3,1,0,0,1.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,1,numNodes,0,0,0.0);
sysEQ.addBC(BC2);

sysEQ.ready()


%% time for the convergence study
rhoMax = 1; rhoMin = 0.2; rhoEss = 0.0;
tEnd = 0.15;
nSteps_exact = 2000;
nSteps_set = [128,64,32,16,8];
whichnode = 3;


gs4 = GS4(rhoMax,rhoMin,rhoEss,1.,1);
gs4.convergence_test(sysEQ,tEnd,nSteps_exact,nSteps_set,whichnode)



