function [ ballisticOut ] = ballisticFlux( eta, xi, Kn )
% calculates Eq (49) from chen2002ballistic for a given location and time

F = (@(mu) mu.*exp(-eta./Kn./mu));
mu_t = 1/Kn.*eta./xi;

    if mu_t >= 0 && mu_t <=1
        ballisticOut = 1/2.*integral(F,mu_t,1);
    else
        ballisticOut = 0;
    end

end