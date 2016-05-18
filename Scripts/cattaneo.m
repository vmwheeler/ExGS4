%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear;
clc;
addpath('../Base','../Elements','../ForceTerms/');

%% Physical and numerical constants
numEle = 40;
numNodes = numEle+1;
Kn = 1.0;
rhoMax = 1.0; rhoMin = 0.8; rhoEss = 0.7;
tEnd = 5.;
numSteps = 10;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,2);
gs4damped = GS4(0.0,0.0,0.0,dt,2);

%% Import (or create) mesh data
% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

%% Initialize vectors
problemIC = zeros(numNodes,1);
problemIC(1) = 1;
xPlot = linspace(0,1,numNodes);
sol = zeros(numNodes,1);
xFlux = zeros(numNodes-1,1);
flux = zeros(numNodes-1,1);

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
    ele(i) = CF_1DL(i,eleNodes, [1,1,Kn^2/3,0], fhandle);
end

%% Set up system of equations
% initiallize
sysEQ = SystemEQ(nodes,ele);


%% set BCs
BC1 = BoundaryCondition(1,1,1,0,0,1.0);
sysEQ.addBC(BC1);
BC2 = BoundaryCondition(2,2,numNodes,0,0,-0.3);
sysEQ.addBC(BC2);


%% Solve!
sysEQ.ready()
for n = 1:numSteps
    fprintf('Timestep #%i \n',n)
    gs4.time_march(sysEQ);
    % find heat flux
    sysEQ.computeDirs()
    for k = 1:numEle
        ele(k).computeFlux(gs4,Kn,0);
        flux(k) = ele(k).qNew;
    end
end

%% Post-process
for j = 1:numNodes
    sol(j) = nodes(j).y;
end

sysEQ.computeDirs()
for i = 1:numEle
    xFlux(i) = (sysEQ.ele(i).nodes(2).loc ...
        + sysEQ.ele(i).nodes(1).loc)/2;
end

sysEQ.reset();

%% Solve!
sysEQ.ready()
for n = 1:numSteps
    fprintf('Timestep #%i \n',n)
    gs4damped.time_march(sysEQ);
    % find heat flux
    sysEQ.computeDirs()
    for k = 1:numEle
        ele(k).computeFlux(gs4damped,Kn,0);
        fluxdamped(k) = ele(k).qNew;
    end
end

%% Post-process
for j = 1:numNodes
    soldamped(j) = nodes(j).y;
end


% Display results
figure(3)
h = plot(xPlot,sol,'-rs',xPlot,soldamped,'--ks');
set(gca,'fontsize',8)
set(gcf, 'Units','centimeters', 'Position',[0 0 6 5])
xData = get(h,'XData');
set(gca,'Xtick',linspace(xData{1}(1),xData{1}(end),5))
ylabel('dimensionless temperature')
xlabel('dimensionless distance')
ylim([0,1])
xlim([0,1])

figure(2)
h = plot(xFlux,flux,'-ob',xFlux,fluxdamped,'--ok');
set(gca,'fontsize',8)
set(gcf, 'Units','centimeters', 'Position',[0 0 6 5])
xData = get(h,'XData');
set(gca,'Xtick',linspace(xData{1}(1),xData{1}(end),5))
ylabel('dimensionless flux')
xlabel('dimensionless distance')

