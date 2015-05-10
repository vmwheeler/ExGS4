% function [ tempOut ] = ballisticTemp( KN, ETA, XI )
% % calculates Eq (50) from chen2002ballistic for a given location and time
% 
% F = (@(mu)exp(-ETA./mu./KN));
% 
% MU_T = (1/KN)*ETA/XI;
% if MU_T >= 0 && MU_T <=1
%     tempOut = 1/2 * ( quadgk(F,MU_T,1) );
%     %tempOut = 1/2 * ( quadgk(F,0,1) );
% else
%    tempOut = 0;
% end
% 
% end

function [ tempOut ] = ballisticTemp( KN, ETA, XI )
% calculates Eq (50) from chen2002ballistic for a given location and time

F = (@(mu) heaviside(XI-ETA./KN./mu).*exp(-ETA./KN./mu));

%MU_T = (1/KN)*ETA/XI;
%if MU_T >= 0 && MU_T <=1
    
   tempOut = 1/2 * ( quadgk(F,0,1,'RelTol',1e-6,'AbsTol',1e-6) );
%else
%   tempOut = 0;
%end

end