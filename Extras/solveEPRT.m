function [ WantedProfile, WantedFlux ] = solveEPRT( ...
    knud, tend, numt, numEl, gamgam, rho, rhoSp )

% gauss points and weights for integration of all angles
%MuVec = [ 1 -1 ];
%MuVec = [ 0.577350269189 -0.577350269189 ];
MuVec = [ 0.095012509837637 -0.095012509837637 0.281603550779258 ...
        -0.281603550779258 0.458016777657227 -0.458016777657227...
        0.617876244402643 -0.617876244402643 0.755404408355003 ...
        -0.755404408355003 0.865631202387831 -0.865631202387831...
        0.944575023073232 -0.944575023073232 0.989400934991649...
        -0.989400934991649];
         
nOrd = length(MuVec);

%Weights = [ 1 1 ];
Weights = [ 0.189450610455068 0.189450610455068 0.182603415044923...  
       0.182603415044923 0.169156519395002 0.169156519395002...
       0.149595988816576 0.149595988816576 0.124628971255533...
       0.124628971255533 0.095158511682492 0.095158511682492...
       0.062253523938647 0.062253523938647 0.027152459411754...
       0.027152459411754];

dist = 1;
nEle = numEl;
h = dist/nEle;
nNodes = nEle+1;

kn = knud;
sigma = 51.4995;
TL = 300.1;
TR = 300;

KfaT = zeros(nNodes*nOrd,nNodes*nOrd);
CfaT = zeros(nNodes*nOrd,nNodes*nOrd);
Ivec = zeros(nNodes*nOrd,1);

% calculate [C] entries (done here for speed)
c11 = h/3;
c12 = h/6;
c21 = h/6;
c22 = h/3;
% note!:  [K] contribution for f (by itself) is the same so use in [K]
% assembly


for j = 1:nEle
    for i = 1:nOrd
        loc = (j-1)*nOrd+i;
        nxt = loc + nOrd;
        mu = MuVec(i);
               
        if mu > 0
            gamma = gamgam;
        else
            gamma = -gamgam;
        end
        
        % assemble [K] matrix
        KfaT(loc,loc) = KfaT(loc,loc) + kn*mu/2*(-1+gamma) + c11;
        KfaT(loc,nxt) = KfaT(loc,nxt) + kn*mu/2*(1-gamma) + c12;
        KfaT(nxt,loc) = KfaT(nxt,loc) + kn*mu/2*(-1-gamma) + c21;
        KfaT(nxt,nxt) = KfaT(nxt,nxt) + kn*mu/2*(1+gamma) + c22;
        
        % assemble [C] matrix
        CfaT(loc,loc) = CfaT(loc,loc) + c11;
        CfaT(loc,nxt) = CfaT(loc,nxt) + c12;
        CfaT(nxt,loc) = CfaT(nxt,loc) + c21;
        CfaT(nxt,nxt) = CfaT(nxt,nxt) + c22;
        
        % integral contribution to [K]
        for k = 1:nOrd
            intCol = (j-1)*nOrd+k;
            intRow = loc;
            intCnxt = intCol + nOrd;
            intRnxt = intRow + nOrd;
            w = Weights(k);
            
            KfaT(intRow,intCol) = KfaT(intRow,intCol) - w*h/6;
            KfaT(intRow,intCnxt) = KfaT(intRow,intCnxt) - w*h/12;
            KfaT(intRnxt,intCol) = KfaT(intRnxt,intCol) - w*h/12;
            KfaT(intRnxt,intCnxt) = KfaT(intRnxt,intCnxt)  - w*h/6;
        end
        %disp(KfaT )
%fprintf('\n **next integration** \n')
%pause
    end
end

% evaluate F
Fnow = zeros(nNodes*nOrd,1);
Fltr = zeros(nNodes*nOrd,1);
Isol = zeros(nNodes*nOrd,1);

% set initial conditions
for i = 1:nNodes*nOrd
    Isol(i) = sigma*(TR)^4/2;
end
for i = 1:nOrd
    mu = MuVec(i);
    if mu > 0
        Isol(i) = sigma*(TL)^4/2;
    end
end
Idot = CfaT \ (-KfaT*Isol);

Itot = zeros(nNodes,1);
Temp = zeros(nNodes,1);

bcLeft = sigma*(TL)^4/2;
bcRght = sigma*(TR)^4/2;

for j = 1:nNodes
    Isum = 0;
    for i = 1:nOrd
        w = Weights(i);
        loc = (j-1)*nOrd+i;
        Isum = Isum + w*Isol(loc);
    end
    Itot(j) = Isum;
    Temp(j) = ((Isum/sigma))^0.25;
end
SolCelL{1,1}=Temp;
SolCelL{2,1}=(Temp-TR)/(TL-TR);

%set parameters for GS4
rho1 = rho;
rho2 = rhoSp;

%set parameters for temporal discretization
endTime = tend;
numSteps = numt;
delt = endTime/numSteps;

for n = 1:numSteps
    
    [ Idot, Isol ] = GSSSSone_EPRT( CfaT, KfaT, Fnow, Fltr, ...
        Idot, Isol, delt, nNodes, rho1, rho2, MuVec, bcLeft, bcRght, nOrd);
    
    %[ Isol ] = GenTrap_EPRT( CfaT, KfaT, Fnow, Fltr, ...
    %Isol, delt, nNodes, 0.0, MuVec, bcLeft, bcRght, nOrd);

    
    %SolCelL{1,n+1}=Temp;
    %SolCelL{2,n+1}=(Temp-TR)/(TL-TR);
    
    t = n * delt;
    
end


q = zeros(nNodes,1);

for j = 1:nNodes
    Isum = 0;
    
    for i = 1:nOrd
        w = Weights(i);
        loc = (j-1)*nOrd+i;
        Isum = Isum + w*Isol(loc);
        q(j) = q(j) + MuVec(i)*w*Isol(loc);
    end
    Itot(j) = Isum;
    Temp(j) = ((Isum/sigma))^0.25;
end


%WantedProfile = SolCelL{2,numSteps};
WantedProfile = (Temp-TR)/(TL-TR);
WantedFlux = q/sigma/(TL^4-TR^4);
end