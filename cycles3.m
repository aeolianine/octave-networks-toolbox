% Calculate number of loops/cycles of length 3
%
% INPUTs: adj - adjacency matrix, nxn
% OUTPUTs: L3 - number of triangles (loops of length 3)
%
% Note: Valid for an undirected network.
% Last updated: Oct 5, 2012

function l3 = cycles3(adj)

    l3 = trace(adj^3) / 6; % trace(adj^3)/3!


%!test
%!shared T
%! T = load_test_graphs();
%!assert(cycles3(T{4}{2}),2)
%!assert(cycles3(T{18}{2}),0)
%!assert(cycles3(T{13}{2}),1)
%!assert(cycles3(edgeL2adj(canonicalNets(randi(10)+3,'btree'))),0)
%!assert(cycles3(edgeL2adj(canonicalNets(4,'trilattice'))),2)


%!demo
%! cycle4 = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]; 
%! cycles3(cycle4)
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! cycles3(bowtie)