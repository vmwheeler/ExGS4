% Starting fresh on solving the BDA and the proposed alterations...

clear all;
close all;
clc;

tic;

%% Physical and numerical constants
numEle = 50;
numNodes = 2*numEle+1;
Kn = 10;
tEnd = 10;
timeSteps = 50;
dt=tEnd/timeSteps;

FT = 0.5;

%get GS4 coefficients
BDAgs4 = GS4(0,0,0,dt,2);
NHEgs4 = GS4(1,0,0,dt,1);
BNHEgs4 = GS4(1,0,0,dt,1);
CFgs4 = GS4(0,0,0,dt,2);
Fgs4 = GS4(1,0,0,dt,1);
Cgs4 = GS4(0,0,0,dt,2);

%% Import (or create) mesh data

% generate nodal locations
nodeLocs = linspace(0,1,numNodes);

% create initial values (might be inconvenient when reading in mesh but 
% we shall see...)
problemIC = zeros(numNodes,1);
%problemIC(1) = 1;
% note the yd initial value is set to zero below

% initialize and number nodes
for i = 1:numNodes
    iv = problemIC(i);
    BDAnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
    NHEnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
    BNHEnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
    CFnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
    Fnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
    Cnodes(i) = Node(i,nodeLocs(i),0,iv,0,0);
end

% assign nodes to elements
for i = 1:numEle
    BDAeleNodes{i} = [BDAnodes(2*i-1);BDAnodes(2*i);BDAnodes(2*i+1)];
    NHEeleNodes{i} = [NHEnodes(2*i-1);NHEnodes(2*i);NHEnodes(2*i+1)];
    BNHEeleNodes{i} = [BNHEnodes(2*i-1);BNHEnodes(2*i);BNHEnodes(2*i+1)];
    CFeleNodes{i} = [CFnodes(2*i-1);CFnodes(2*i);CFnodes(2*i+1)];
    FeleNodes{i} = [Fnodes(2*i-1);Fnodes(2*i);Fnodes(2*i+1)];
    CeleNodes{i} = [Cnodes(2*i-1);Cnodes(2*i);Cnodes(2*i+1)];
end

%create elements
for i = 1:numEle
    BDAele(i) = BDAeleQ(i,BDAeleNodes{i}, [1,1,Kn^2/3],[0,0,0],BDAgs4);
    NHEele(i) = Element1DQ(i,NHEeleNodes{i}, [1,0],[0,0,0],NHEgs4);
    BNHEele(i) = Element1DQ(i,BNHEeleNodes{i}, [1,0],[0,0,0],BNHEgs4);
    CFele(i) = CFeleQ(i,CFeleNodes{i}, [1,1,Kn^2/3,FT],[0,0,0],CFgs4);
    Fele(i) = Element1DQ(i,FeleNodes{i}, [1,Kn^2/3],[0,0,0],Fgs4);
    Cele(i) = CFeleQ(i,CeleNodes{i}, [1,1,Kn^2/3,0],[0,0,0],Cgs4);
end


%% Construct Global system

% initiallize
BDAsys = SystemEQ(numNodes,numEle,BDAnodes,BDAgs4);
NHEsys = SystemEQ(numNodes,numEle,NHEnodes,NHEgs4);
BNHEsys = SystemEQ(numNodes,numEle,BNHEnodes,BNHEgs4);
CFsys = SystemEQ(numNodes,numEle,CFnodes,CFgs4);
Fsys = SystemEQ(numNodes,numEle,Fnodes,Fgs4);
Csys = SystemEQ(numNodes,numEle,Cnodes,Cgs4);

% add local stiffness matrices to global system
for i = 1:numEle
   BDAsys.addElement(BDAele(i));
   NHEsys.addElement(NHEele(i));
   BNHEsys.addElement(BNHEele(i));
   CFsys.addElement(CFele(i));
   Fsys.addElement(Fele(i));
   Csys.addElement(Cele(i));
end

%%BC1
BC1bda = BoundaryCondition(1,4,1,0.5*Kn,0.5*Kn,0);
BDAsys.addBC(BC1bda)

BC1bnhe = BoundaryCondition(1,3,1,0,0.5*Kn,0);
BNHEsys.addBC(BC1bnhe)

BC1 = BoundaryCondition(1,3,1,0,0.5*Kn,0.5*Kn);
NHEsys.addBC(BC1)

BC1cf = BoundaryCondition(1,4,1,0.5*Kn,0.5*Kn,0.5*Kn);
CFsys.addBC(BC1cf)

BC1f = BoundaryCondition(1,3,1,0,0.5*Kn,0.5*Kn);
%BC1f = BoundaryCondition(1,1,1,0,0,1);
Fsys.addBC(BC1f)

BC1c = BoundaryCondition(1,4,1,0.5*Kn,0.5*Kn,0.5*Kn);
%BC1c = BoundaryCondition(1,1,1,0,0,1);
Csys.addBC(BC1c)


%% BC2
BC2bda = BoundaryCondition(2,4,numNodes,0.5*Kn,0.5*Kn,0);
BDAsys.addBC(BC2bda)

BC2bnhe = BoundaryCondition(2,3,numNodes,0,0.5*Kn,0);
BNHEsys.addBC(BC2bnhe)

BC2 = BoundaryCondition(2,3,numNodes,0,0.5*Kn,0);
NHEsys.addBC(BC2)

BC2cf = BoundaryCondition(2,4,numNodes,0.5*Kn,0.5*Kn,0);
CFsys.addBC(BC2cf)

BC2f = BoundaryCondition(2,3,numNodes,0,0.5*Kn,0);
%BC2f = BoundaryCondition(2,1,numNodes,0,0,0);
Fsys.addBC(BC2f)

BC2c = BoundaryCondition(2,4,numNodes,0.5*Kn,0.5*Kn,0);
%BC2c = BoundaryCondition(2,1,numNodes,0,0,0);
Csys.addBC(BC2c)

BDAsys.setBCs()
NHEsys.setBCs()
BNHEsys.setBCs()
CFsys.setBCs()
Fsys.setBCs()
Csys.setBCs()

% initialize garbage for qBDA calc
qOld = zeros(numEle,1);
qNew = zeros(numEle,1);
qdOld = zeros(numEle,1);
qdNew = zeros(numEle,1);
ypOld = zeros(numEle,1);

qcOld = zeros(numEle,1);
qcNew = zeros(numEle,1);
qcdOld = zeros(numEle,1);
qcdNew = zeros(numEle,1);
ypcOld = zeros(numEle,1);

qcfOld = zeros(numEle,1);
qcfNew = zeros(numEle,1);
qcfdOld = zeros(numEle,1);
qcfdNew = zeros(numEle,1);
ypcfOld = zeros(numEle,1);
ypcfdOld = zeros(numEle,1);


% gs4 parameters for q calc
W1 = BDAgs4.w1;
L5W2 = BDAgs4.lam5w2;
L6W1 = BDAgs4.lam6w1;
l5 = BDAgs4.lam5;
delt = BDAgs4.dt;

t=0;
for j = 1:timeSteps  
    
    t = t + dt;
    tw = t + BDAgs4.w1*dt;
    
    
    %add ballistic force
    for q = 1:numEle
        x1 = BDAsys.ele{q}.nodes{1}.loc;
        x2 = BDAsys.ele{q}.nodes{2}.loc;
        x3 = BDAsys.ele{q}.nodes{3}.loc;
        
        BDAsys.ele{q}.force(1) = ballisticTerm(Kn,x1,tw);
        BDAsys.ele{q}.force(2) = ballisticTerm(Kn,x2,tw);
        BDAsys.ele{q}.force(3) = ballisticTerm(Kn,x3,tw);
                
        NHEsys.ele{q}.const(2) = (1-exp(-tw))*Kn^2/3;
        
        BNHEsys.ele{q}.const(2) = (1-exp(-tw))*Kn^2/3;
        
        BNHEsys.ele{q}.force(1) = ballisticTemp(Kn,x1,tw);
        BNHEsys.ele{q}.force(2) = ballisticTemp(Kn,x2,tw);
        BNHEsys.ele{q}.force(3) = ballisticTemp(Kn,x3,tw);
    end
    
    BDAsys.updateSystem()
    NHEsys.updateSystem()
    BNHEsys.updateSystem()
    
    BDAsys.solve(10E-7,10);
    BDAsys.timeMarch();
    
    % find heat flux for hyperbolic models (messy)
    BDAsys.computeDirs()
    Csys.computeDirs()
    CFsys.computeDirs()
    for k = 1:numEle
        qdNew(k) = ( -Kn/3*(ypOld(k)+W1*(BDAsys.ele{k}.yp-ypOld(k)))...
                     +(L6W1+L5W2*delt-W1*delt-1)*qdOld(k) ...
                     - qOld(k))/(L6W1+L5W2*delt);
        
        qNew(k) = qOld(k) + delt*qdOld(k) ...
            + l5*delt*(qdNew(k)-qdOld(k)); 
        
        qcdNew(k) = ( -Kn/3*(ypcOld(k)+W1*(Csys.ele{k}.yp-ypcOld(k)))...
                     +(L6W1+L5W2*delt-W1*delt-1)*qcdOld(k) ...
                     - qcOld(k))/(L6W1+L5W2*delt);
        
        qcNew(k) = qcOld(k) + delt*qcdOld(k) ...
            + l5*delt*(qcdNew(k)-qcdOld(k));
        
        
        qcfdNew(k) = ( -Kn/3*(ypcfOld(k)+W1*(CFsys.ele{k}.yp-ypcfOld(k)))...
                       -Kn/3*FT*(ypcfdOld(k)+W1*(CFsys.ele{k}.ydp-ypcfdOld(k)))...
                     +(L6W1+L5W2*delt-W1*delt-1)*qcfdOld(k) ...
                     - qcfOld(k))/(L6W1+L5W2*delt);
        
        qcfNew(k) = qcfOld(k) + delt*qcfdOld(k) ...
            + l5*delt*(qcfdNew(k)-qcfdOld(k)); 
    end
    
    NHEsys.solve(10E-7,10);
    NHEsys.timeMarch();
    
    BNHEsys.solve(10E-7,10);
    BNHEsys.timeMarch();
    
    CFsys.solve(10E-7,10);
    CFsys.timeMarch();
    
    Fsys.solve(10E-7,10);
    Fsys.timeMarch();
    
    Csys.solve(10E-7,10);
    Csys.timeMarch();
    
    qOld = qNew;
    qdOld = qdNew;
    qcOld = qcNew;
    qcdOld = qcdNew;
    qcfOld = qcfNew;
    qcfdOld = qcfdNew;
    for m = 1:numEle
        ypOld(m) = BDAsys.ele{m}.yp;
        ypcOld(m) = Csys.ele{m}.yp;
        ypcfOld(m) = CFsys.ele{m}.yp;
        ypcfdOld(m) = CFsys.ele{m}.ydp;
    end
end



%% Post-process
    
BDAsys.computeDirs()
NHEsys.computeDirs()
BNHEsys.computeDirs()
CFsys.computeDirs()
Fsys.computeDirs()
Csys.computeDirs()

for i = 1:numEle
    %xFlux(i) = (BDAsys.ele{i}.nodes{1}.loc+BDAsys.ele{i}.nodes{2}.loc)/2;
    xFlux(i) = BDAsys.ele{i}.nodes{2}.loc;
    qNHE(i) = -(1-exp(-t))*Kn/3*NHEsys.ele{i}.yp;
    qF(i) = -Kn/3*Fsys.ele{i}.yp;
    BNHEflux(i) = -(1-exp(-t))*Kn/3*BNHEsys.ele{i}.yp;
    ballFlux(i) = ballisticFlux(Kn, xFlux(i), t);
end

qBDA = qNew' + ballFlux;
qBNHE = BNHEflux + ballFlux;

% solve EPRT
eprtEle = 100;
rhoMinEPRT = 0.0;
rhoEssEPRT = 0.0;
gamgam = 1;
[Teprt, qEPRT] = solveEPRT( Kn, tEnd, round(50), eprtEle, gamgam,...
    rhoMinEPRT, rhoEssEPRT );

xEPRT = linspace(0,1,eprtEle+1);


scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/5 scrsz(4)])

positionVector1 = [0.15, 0.03, 0.8, 0.44];
subplot(2,1,2,'Position',positionVector1)
qs = plot(xEPRT, qEPRT,...
    xFlux, qBDA, ...
    xFlux, qNHE,...
    xFlux, qcfNew,...
    xFlux, qF,...
    xFlux, qcNew);
ylabel('Nondimensional Flux')
xlabel('Nondimensional Coordinate')

%axis([0 1 -0.1 0.3])


BDAdiff = zeros(numNodes,1);
BNHEdiff = zeros(numNodes,1);
for j = 1:numNodes
    BDAdiff(j) = BDAnodes(j).yNew;
    BNHEdiff(j) = BNHEnodes(j).yNew;
    Tnhe(j) = NHEnodes(j).yNew;
    Tcf(j) = CFnodes(j).yNew;
    Tf(j) = Fnodes(j).yNew;
    Tc(j) = Cnodes(j).yNew;
end

xPlot = linspace(0,1,numNodes);

Tball = zeros(numNodes,1);
for j = 1:numNodes
    Tball(j) = ballisticTemp( Kn, xPlot(j), t );
end

Tbda = BDAdiff + Tball;
Tbnhe = BNHEdiff + Tball;

positionVector1 = [0.15, 0.5, 0.8, 0.44];
subplot(2,1,1,'Position',positionVector1)
curves = plot(xEPRT,Teprt,...
    xPlot,Tbda,...
    xPlot,Tnhe,...
    xPlot,Tcf,...
    xPlot,Tf,...
    xPlot,Tc);
title(strcat('t*=',num2str(t),',   ','Kn=',num2str(Kn)))
ylabel('Nondimensional Temperature')
axis([0 1 0 1])
legend('EPRT','BDA','NHE','C-F','Fourier','Cattaneo')


%cd fig_out/
%write solution to file for nice plotting
eprtMat = [xEPRT' , Teprt ];
restMat = [xPlot', Tbda, Tnhe',Tcf',Tf',Tc'];
filename = strcat('fig_out2/eprt','Kn',num2str(100*Kn),'t',num2str(100*tEnd),'.dat');
filename2 = strcat('fig_out2/','Kn',num2str(100*Kn),'t',num2str(100*tEnd),'.dat');
dlmwrite(filename,eprtMat);
dlmwrite(filename2,restMat);


eprtMatQ = [xEPRT' , qEPRT ];
restMatQ = [xFlux',qBDA',qNHE',qcfNew,qF',qcNew];
filenameQ = strcat('fig_out2/Qeprt','Kn',num2str(100*Kn),'t',num2str(100*tEnd),'.dat');
filenameQ2 = strcat('fig_out2/Q','Kn',num2str(100*Kn),'t',num2str(100*tEnd),'.dat');
dlmwrite(filenameQ,eprtMatQ);
dlmwrite(filenameQ2,restMatQ);


fprintf('computation time = %f\n', toc)
