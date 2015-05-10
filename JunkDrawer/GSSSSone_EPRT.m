function [ VoneNew, VtwoNew ] = GSSSSone_EPRT( MonE, MtwO, Fnow, Fltr, ...
    Vone, Vtwo, delt, nNodes, rho1, rho2, MuVec, bcLeft, bcRght, nOrd)


lam6W1 = (3+rho1+rho2-rho1*rho2)/(2*(1+rho1)*(1+rho2));
lam5W2 = 1/((1+rho1)*(1+rho2));
lam4W1 = 1/(1+rho1);
w1 = 1/(1+rho1);
lam4 = 1;
lam5 = 1/(1+rho2);

LeftMaT = lam6W1*MonE + lam5W2*delt*MtwO;

%disp(LeftMaT)

% apply boundary conditions
Rhs = -MonE*Vone - MtwO*(Vtwo + lam4W1*delt*Vone)...
    + (1 - w1)*Fnow + w1*Fltr;


for i = 1:nOrd
    mu = MuVec(i);
    
    if mu > 0
        LeftMaT(i,:) = zeros;
        LeftMaT(i,i) = 1;
        DelvBnd = (bcLeft - (Vtwo(i)+ lam4*Vone(i)*delt))/(lam5*delt);
        Rhs(i) = DelvBnd;
    end
end

j = nNodes;
for i = 1:nOrd
    loc = (j-1)*nOrd+i;
    mu = MuVec(i);
    
    if mu < 0
        LeftMaT(loc,:) = zeros;
        LeftMaT(loc,loc) = 1;
        DelvBnd = (bcRght - (Vtwo(loc)+ lam4*Vone(loc)*delt))/(lam5*delt);
        Rhs(loc) = DelvBnd;
    end
end


%disp(Rhand);

Delv = LeftMaT \ Rhs;

VoneNew = Vone + Delv;
VtwoNew = Vtwo + lam4*Vone*delt + lam5*delt*Delv;



end
