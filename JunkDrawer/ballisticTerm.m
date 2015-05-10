% function [ ballisticOut ] = ballisticTerm( KN, ETA, XI )
% % calculates Eq (50) from chen2002ballistic for a given location and time
% 
% F = (@(mu) exp(-ETA./mu./KN));
% 
% MU_T = 1/KN*ETA/XI;
% if MU_T >= 0 && MU_T <=1
%    ballisticOut = 1/2/KN * ( quadgk(F,MU_T,1) + ETA/(KN*XI^2)*exp(-XI) );
%    %ballisticOut = 1/2/KN * ( quadgk(F,MU_T,1) );
% else
%    ballisticOut = 0;
% end
% 
% 
% end

% function [ ballisticOut ] = ballisticTerm( KN, ETA, XI )
% % calculates Eq (50) from chen2002ballistic for a given location and time
% 
% F = (@(mu) heaviside(XI-ETA./KN./mu).*exp(-ETA./KN./mu));
% 
% MU_T = 1/KN*ETA/XI;
% if MU_T >= 0 && MU_T <=1
%    ballisticOut = 1/2/KN * ( quadgk(F,0,1) + (ETA^2./KN./XI.^2).*exp(-XI) );
% else
%    ballisticOut = 0;
% end
% 
% 
% end

function [ ballisticOut ] = ballisticTerm( KN, ETA, XI )
% calculates Eq (50) from chen2002ballistic for a given location and time

%tic;

F = (@(mu) heaviside(XI-ETA./KN./mu).*exp(-ETA./KN./mu));

MU_T = 1/KN*ETA/XI;

ballisticOut = 1/2/KN *  quadgk(F,0,1,'RelTol',1e-6,'AbsTol',1e-6);

if MU_T >= 0 && MU_T <=1
   ballisticOut = ballisticOut + 1/2 * (ETA./KN^2./XI.^2).*exp(-XI) ; 
end

ballisticOut = KN * ballisticOut;

%fprintf('int time = %f\n',toc)

end