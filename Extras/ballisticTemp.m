function [ tempOut ] = ballisticTemp( eta, xi, Kn )
% calculates Eq (51) from chen2002ballistic for a given location and time

    F = (@(mu) exp(-eta./Kn./mu));

    mu_t = (1/Kn)*eta/xi;
    if mu_t >= 0 && mu_t <=1    
       tempOut = 1/2 * ( integral(F,mu_t,1) );
    else
       tempOut = 0;
    end

end