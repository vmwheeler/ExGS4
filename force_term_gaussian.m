function [ output ] = force_term_gaussian( x )
% This is the forcing function that will be integrated in the element
% object and applied to the nodal equations
a = 0.01;
b = 0.5;
c = 0.015;
output = a*exp(-(x - b).^2./c);

end

