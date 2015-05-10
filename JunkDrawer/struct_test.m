element1 = struct('node_locs', [0.2,0.35], 'ele_num', 3, ...
    'stiffness', [1,-1;-1,1]);
element2 = struct('node_locs', [0.2,0.35], 'ele_num', 3, ...
    'stiffness', [1,-1;-1,1]);
elements = [element1,element2,'teststring'];

systemEQ = struct('elements')

%okay!  time to reformat the who shebang.