function [ element ] = EleCreateLinear1D( num, node_in, const_in, f_in)
% This function creates a 1D linear element

dx = abs(node_in(1).location - node_in(2).location);

kmat = 1/dx * [ 1 -1 ; ...
               -1  1 ];
               
cmat = dx/6 * [ 2 1 ; ...
                1 2 ];   
                
mmat = cmat;

element = struct('num', num, ...
                 'nodes', node_in, ...
                 'constant', const_in, ...
                 'force', f_in, ...
                 'm_matrix', mmat, ...
                 'c_matrix', cmat, ...
                 'k_matrix', kmat);


end
