%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear all;
close all;
clc;

%% Physical and numerical constants
numEle = 40;
numNodes = numEle+1;
Kn = 0.1;
rhoMax = 1;
rhoMin = 0.;
rhoEss = 0.;
tEnd = 10000;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

%% Initialize vectors
problemIC = zeros(numNodes,1);
%ele(numEle) = Element1DL();
sol = zeros(numNodes,1);
flux = zeros(numNodes-1,1);

%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

problemIC(1) = 1;
% initialize and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
end

%fhandle = @force_term_0;
fhandle = @force_term_gaussian;

% assign nodes to elements
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = Element1DL(i,eleNodes, [1,Kn^2/3], fhandle);
end

% initiallize
sysEQ = SystemEQ(numNodes,numEle,nodes,gs4,ele(1));

% construct system
for q = 1:numEle
    sysEQ.addElement(ele(q));
end

% set BCs
BC1 = BoundaryCondition(1,1,1,0,0,0.2);
%BC1 = BoundaryCondition(1,3,1,0,0.5*Kn,0.5*Kn);
%BC1 = BoundaryCondition(1,2,1,0,0,0.001);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,1,numNodes,0,0,0.0);
%BC2 = BoundaryCondition(2,3,numNodes,0,0.5*Kn,0);
%BC2 = BoundaryCondition(2,2,numNodes,0,0,-0.001);
sysEQ.addBC(BC2);


for n = 1:numSteps
    fprintf('Timestep #%i \n',n)
    gs4.time_march(sysEQ);
end

xPlot = linspace(0,1,numNodes);
for j = 1:numNodes
    sol(j) = nodes(j).yNew;
end

sysEQ.computeDirs()
for i = 1:numEle
    xFlux(i) = sysEQ.ele(i).nodes(2).loc;
    flux(i) = -Kn^2/3*sysEQ.ele(i).yp;
end

figure(3)
plot(xPlot,sol,'-rs',xPlot,-sysEQ.force,'--b')
ylabel('u')
%ylim([0,1])
%xlim([0,1])

figure(2)
plot(xFlux,flux,'-ob')
ylabel('flux')

