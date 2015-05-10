function [ node ] = NodeCreate( num, loc, on_or_off )
% This function creates a node by making a structure array that holds all
% relevant nodal information

node = struct('num',num, ...
              'location', loc, ...
              'y_new', 0, ...
              'yd_new', 0, ...
              'ydd_new', 0, ...
              'y_old', 0, ...
              'yd_old', 0, ...
              'ydd_old', 0, ...
              'on_boundary', on_or_off);


end

