function [ output ] = force_term_gaussian( x, t, a, b, c )
% This is the forcing function that will be integrated in the element
% object and applied to the nodal equations

output = a*exp(-(x - b).^2./c);

end

