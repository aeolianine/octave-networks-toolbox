% The longest shortest path between any two nodes nodes in the network.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: network diameter
%
% Other routines used: simpleDijkstra.m
% Last updated: Oct 8 2012

function diam = diameter(adj)

    diam = 0;

    for i = 1:size(adj, 1)
        d = simpleDijkstra(adj, i);
        diam = max([max(d), diam]);
    end


%!test
%!shared T, el, adj
%! T = load_test_graphs();
%!assert(diameter(T{13}{2}),1)
%!assert(diameter(T{4}{2}),3)

%! el=canonicalNets(randi(10)+5,'line');
%! adj = edgeL2adj(el);
%!assert(diameter(adj),length(adj)-1)

%! el=canonicalNets(randi(10)+5,'cycle');
%! adj = edgeL2adj(el);
%!assert(diameter(adj),floor(length(adj)/2))

%!demo
%! diameter( [0 1 1; 1 0 1; 1 1 0] )  % cycle-3
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! diameter(bowtie)
%! diameter([0 1; 0 0])