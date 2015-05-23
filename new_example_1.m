%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear all;
close all;
clc;

%% Physical and numerical constants
numEle = 3;
numNodes = numEle+1;
Kn = 0.1;
rhoMax = 1;
rhoMin = 0;
rhoEss = 0;
tEnd = 22;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

% create initial values (might be inconvenient when reading in mesh but 
% we shall see...)
problemIC = zeros(numNodes,1);
problemIC(1) = 1;
% note the yd initial value is set to zero below

% initialize and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    %nodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
    nodes(i) = Node(i,nodeLocs(i),0,0,iv,0);
end

% assign nodes to elements
%ele(numEle) = Element1DL();
for i = 1:numEle
    eleNodes = [nodes(i);nodes(i+1)];
    ele(i) = Element1DL(i,eleNodes, [1,Kn^2/3],[0,0,0],gs4);
end


% initiallize
sysEQ = SystemEQ(numNodes,numEle,nodes,gs4,ele(1));

% construct system
for q = 1:numEle
    sysEQ.addElement(ele(q));
end

% set BCs
BC1 = BoundaryCondition(1,1,1,0,0,1);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,1,numNodes,0,0,0);
sysEQ.addBC(BC2);


for n = 1:numSteps
    
    gs4.time_march(sysEQ);
    
end

xPlot = linspace(0,1,numNodes);
for j = 1:numNodes
    sol(j) = nodes(j).ydNew;
    %sol(j) = nodes(j).yNew;
end

sysEQ.computeDirs()
for i = 1:numEle
    xFlux(i) = sysEQ.ele(i).nodes(2).loc;
    flux(i) = -Kn^2/3*sysEQ.ele(i).yp;
end

figure(3)
plot(xPlot,sol,'-rs')
ylabel('u')

figure(2)
plot(xFlux,flux,'-ob')
ylabel('flux')

