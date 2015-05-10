function [ ballisticOut ] = ballisticFlux( KN, ETA, XI )
% calculates Eq (50) from chen2002ballistic for a given location and time

%tic;

F = (@(mu) mu.*heaviside(XI-ETA./KN./mu).*exp(-ETA./KN./mu));

ballisticOut = 1/2 *  quadgk(F,0,1,'RelTol',1e-6,'AbsTol',1e-6);

ballisticOut = ballisticOut;

%fprintf('int time = %f\n',toc)

end