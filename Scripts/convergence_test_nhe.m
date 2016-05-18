%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear;
close all;
clc;
addpath('../Base','../Elements','../ForceTerms/');

%% Physical and numerical constants
numEle = 10;
numNodes = numEle+1;
Kn = 1.0;


%% Initialize vectors
problemIC = zeros(numNodes,1);
sol = zeros(numNodes,1);
flux = zeros(numNodes-1,1);
xPlot = linspace(0,1,numNodes);

%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

problemIC(1) = 1;
% initialize and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end

%% Create elements from mesh
% set the force term (simply 0 here)
fhandle = @force_term_0;

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = CF_1DL(i,eleNodes, [1,1,Kn^2/3,1], fhandle);
end

%% Set up system of equations
% initiallize
sysEQ = SystemEQ(nodes,ele);

%% set BCs
BC1 = BoundaryCondition(1,3,1,0,Kn*0.5,Kn*0.5);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0,Kn*0.5,0.0);
sysEQ.addBC(BC2);

%% Solve!
% grab the original bigK so we can adjust the time-varying coefficient
bigK0 = sysEQ.bigK;
sysEQ.ready()

%% time for the convergence study
rhoMax = 1; rhoMin = 0.44; rhoEss = 0.32;
tEnd = 0.15;
nSteps_exact = 1000;
nSteps_set = [128,64,32,16,8];
whichnode = 3;

gs4 = GS4(rhoMax,rhoMin,rhoEss,1.,1);
gs4.convergence_test(sysEQ,tEnd,nSteps_exact,nSteps_set,whichnode)



