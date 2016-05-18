T%A simple example of using ExGS4 to solve a 1D heat conduction problem
clear;
clc;
addpath('../Base','../Elements','../ForceTerms/');

%% Physical and numerical constants
numEle = 20;
numNodes = numEle+1;
Kn = 0.5;
gamma = 0.8;
FTs = [0.0, 0.5, 1.0];
rhoMax = 1.0; rhoMin = 0.2; rhoEss = 0.2;
tEnd = 2.0;
numSteps = 20;
dt=tEnd/numSteps;

gs4 = GS4(rhoMax,rhoMin,rhoEss,dt,2);

%% Import (or create) mesh data
% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

%% Initialize vectors
problemIC = zeros(numNodes,1);
xPlot = linspace(0,1,numNodes);
sol = zeros(numNodes,length(FTs));
xFlux = zeros(numNodes-1,1);
flux = zeros(numNodes-1,length(FTs));

% initialize and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    nodes(i) = Node(i,nodeLocs(i),iv,0);
end

% set the force term (simply 0 here)
fhandle = @force_term_0;

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
ylim([-0.1,2.5])
xlim([0,1])

for j = 1:length(FTs)
    FT = FTs(j);

    % assign nodes to elements
    for i = 1:numEle
        eleNodes = [nodes(i);nodes(i+1)];
        ele(i) = CF_1DL(i,eleNodes, [1,1,Kn^2/3,FT], fhandle);
    end

    % initiallize
    sysEQ = SystemEQ(nodes,ele);

    % set BCs
    %% set BCs
    BC1 = BoundaryCondition(1,3,1,Kn*gamma,Kn*gamma,Kn*gamma);
    sysEQ.addBC(BC1);
    BC2 = BoundaryCondition(2,3,numNodes,Kn*gamma,Kn*gamma,0.0);
    sysEQ.addBC(BC2);

    %Time march!
    sysEQ.ready()
    for n = 1:numSteps
        fprintf('Timestep #%i \n',n)
        gs4.time_march(sysEQ);
        sysEQ.computeDirs()
        % find heat flux
        for k = 1:numEle
            ele(k).computeFlux(gs4,Kn,FT);
            flux(k,j) = ele(k).qNew;
        end
    end

    for p = 1:numNodes
        sol(p,j) = nodes(p).y;
    end

    sysEQ.computeDirs()
    for i = 1:numEle
        xFlux(i) = (sysEQ.ele(i).nodes(2).loc ...
            + sysEQ.ele(i).nodes(1).loc)/2;
    end

    sysEQ.reset();
    
    set(0,'CurrentFigure',hT)
    hold on;
    plot(xPlot,sol(:,j), '-s', 'DisplayName',strcat('FT=',num2str(FT,'%2.1f')))
    
    
    set(0,'CurrentFigure',hQ)
    hold on;
    plot(xFlux,flux(:,j), '-o', 'DisplayName',strcat('FT=',num2str(FT,'%2.1f')))
    
end
set(0,'CurrentFigure',hT)
legend(gca,'show')
set(0,'CurrentFigure',hQ)
%legend(gca,'show')
