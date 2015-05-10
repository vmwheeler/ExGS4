% Testing the use of quadratic elements in WOOFE
% CHECK

clear all;
%close all;
clc;

tic;

%% Physical and numerical constants
numEle = 10;
numNodes = 2*numEle+1;
Kn = 0.1;
rhoMax = 1;
rhoMin = 0;
rhoEss = 0;
tEnd = 5;
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
    qnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
end

% assign nodes to elements
for i = 1:numEle
    QeleNodes{i} = [qnodes(2*i-1);qnodes(2*i);qnodes(2*i+1)];
end

%create elements
for i = 1:numEle
    ele(i) = Element1DQ(i,QeleNodes{i}, [1,Kn^2/3],[0,0,0],gs4);
end

% initiallize
sysEQ = SystemEQ(numNodes,numEle,qnodes,gs4);

% construct system
for q = 1:numEle
    sysEQ.addElement(ele(q));
end

% set BCs
BC1 = BoundaryCondition(1,1,1,0,0,1);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,1,numNodes,0,0,0);
sysEQ.addBC(BC2);

sysEQ.setBCs();


for n = 1:numSteps
    
    sysEQ.solve(10E-8,10);
    sysEQ.timeMarch();
    
end

xPlot = linspace(0,1,numNodes);
for j = 1:numNodes
    sol(j) = qnodes(j).yNew;
end

figure(3)
plot(xPlot,sol)

% sysEQ.computeDirs()
% for i = 1:numEle
%     xFlux(i) = sysEQ.ele{i}.nodes{2}.loc;
%     flux(i) = -Kn^2/3*sysEQ.ele{i}.yp;
% end
% 
% figure(2)
% plot(xFlux,flux)
