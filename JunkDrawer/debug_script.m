fprintf('here we go\n')

loc = 0.5;
y = 0.29;
node1 = NodeCreate(1,loc,true)
node2 = NodeCreate(2,0.25,false)
nodes = [node1,node2];

const = [0,1,2];
f = [0,0];

test_element = EleCreateLinear1D(1,nodes,const,f)
