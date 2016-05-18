function [ ballisticOut ] = ballisticTerm( eta, xi, Kn )
    % calculates Eq (50) from chen2002ballistic for a given location and time

    F = @(mu) exp(-eta./Kn./mu);
    mu_t = 1/Kn.*eta./xi;

    if mu_t >= 0 && mu_t <=1
        ballisticOut = 1/2 .*  (integral(F,mu_t,1) +(eta./Kn./xi.^2).*exp(-xi));
    else
        ballisticOut = 0;
    end

end