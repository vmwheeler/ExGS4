%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear all;
close all;
clc;
addpath('../Base','../Elements','../ForceTerms/');

%% Physical and numerical constants
numEle = 10;
numNodes = numEle+1;
Kn = 1.0;
rhoMax = 1; rhoMin = 0.0; rhoEss = 0.0;
tEnd = 1.0;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

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
    nodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
end

% set the force term (simply 0 here)
fhandle = @force_term_0;
%fhandle = @force_term_gaussian;

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
BC1 = BoundaryCondition(1,3,1,0,Kn*0.5,Kn*0.5);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,3,numNodes,0,Kn*0.5,0.0);
sysEQ.addBC(BC2);

bigK0 = sysEQ.bigK;
gs4.tnpw1 = gs4.w1*dt;

fprintf('dt = %f\n', gs4.dt)
for n = 1:numSteps
    fprintf('Timestep #%i \n',n)
    fprintf('gs4.n = %i\n', gs4.n)
    fprintf('gs4.tn = %f\n', gs4.tn)
    fprintf('gs4.tnp1 = %f\n', gs4.tnp1)
    fprintf('gs4.tnpw1 = %f\n', gs4.tnpw1);
    sysEQ.bigK = (1-exp(-gs4.tnpw1))*bigK0;
    gs4.time_march(sysEQ);
end

for j = 1:numNodes
    sol(j) = nodes(j).yNew;
end

sysEQ.computeDirs()
for i = 1:numEle
    xFlux(i) = sysEQ.ele(i).nodes(2).loc;
    flux(i) = -(1-exp(-gs4.tnpw1))*Kn^2/3*sysEQ.ele(i).yp;
end

figure(3)
plot(xPlot,sol,'-rs',xPlot,-sysEQ.force,'--b')
ylabel('dimensionless temperature')
xlabel('dimensionless distance')
ylim([0,1])
xlim([0,1])

figure(2)
plot(xFlux,flux,'-ob')
ylabel('dimensionless flux')
xlabel('dimensionless distance')

