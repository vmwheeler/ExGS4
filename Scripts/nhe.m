%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear;
clc;
addpath('../Base','../Elements','../ForceTerms/');

%% Physical and numerical constants
numEle = 25;
numNodes = numEle+1;
Kn = 0.5;
rhoMax = 1.0; rhoMin = 0.0; rhoEss = 0.0;
tEnd = 5.0;
numSteps = 50;
dt=tEnd/numSteps;
ns = [10,25,50];
tids = ns*dt;


gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,1);

%% Initialize vectors
problemIC = zeros(numNodes,1);
xPlot = linspace(0,1,numNodes);
sol = zeros(numNodes,length(ns));
xFlux = zeros(numNodes-1,1);
flux = zeros(numNodes-1,1);
fluxBD = zeros(numNodes-1,length(ns));

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
    ele(i) = Element1DL(i,eleNodes, [1,Kn^2/3], fhandle);
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
for n = 1:numSteps
    fprintf('Timestep #%i \n',n)
    sysEQ.bigK = (1-exp(-gs4.tnpw1))*bigK0;
    gs4.time_march(sysEQ);
    
    for m = 1:length(ns)
        thisn = ns(m);
        if thisn == n
            for j = 1:numNodes
                sol(j,m) = nodes(j).y;
            end
            sysEQ.computeDirs()
            for i = 1:numEle
                flux(i,m) =  -(1-exp(-gs4.tnpw1))*Kn/3*sysEQ.ele(i).yp;
            end
        end
    end
end

sysEQ.computeDirs()
for i = 1:numEle
    xFlux(i) = (sysEQ.ele(i).nodes(2).loc + sysEQ.ele(i).nodes(1).loc)/2;
end

hT = figure(3);
set(gca,'fontsize',8)
set(gcf, 'Units','centimeters', 'Position',[0 0 6 5])
ylabel('dimensionless temperature')
xlabel('dimensionless distance')
ylim([0,1])
xlim([0,1])
    
hQ = figure(2);
set(gca,'fontsize',8)
set(gcf, 'Units','centimeters', 'Position',[0 0 6 5])
ylabel('dimensionless flux')
xlabel('dimensionless distance')
xlim([0,1])

%% Get EPRT solutions for comparison and plot each snapshot
for i = 1:length(ns)
    [eprtSol,eprtFlux] = solveEPRT(Kn, tids(i), numSteps, numEle, 1.0, 0.0, 0.0 );
    set(0,'CurrentFigure',hT)
    hold on;
    l(i)=plot(xPlot,sol(:,i), '-s', 'DisplayName',strcat('t=',num2str(tids(i),'%3.2f')));
    le = plot(xPlot,eprtSol, '-k', 'DisplayName','EPRT');
    
    set(0,'CurrentFigure',hQ)
    hold on;
    plot(xFlux,flux(:,i), '-o', 'DisplayName',strcat('t=',num2str(tids(i),'%3.2f')))
    plot(xPlot,eprtFlux, '-k')
end
l(end+1) = le;
set(0,'CurrentFigure',hT)
legend(l)
set(0,'CurrentFigure',hQ)

