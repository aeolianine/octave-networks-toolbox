% Calculate the number of independent cycles (use L = m-n+c)
%      where L = num cycles, m - num edges, n - num nodes,
%                      c - number of connected components
% This is also known as the "cyclomatic number": the number of edges
% that needs to be removed so that the graph doesn't have cycles.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: number of independent cycles (or cyclomatic number)
%
% Other routines used: numNodes.m, numEdges.m, findConnComp.m
% Last updated: Oct 5 2012

function L = numCycles(adj)

    n = numNodes(adj);
    m = numEdges(adj);
    comp_mat = findConnComp(adj);

    L = m - n + length(comp_mat);


%!test
%!shared T
%! T = load_test_graphs();
%!assert(numCycles(T{13}{2}),1)
%!assert(numCycles(T{4}{2}),2)
%!assert(numCycles(edgeL2adj(T{10}{2})),0)
%!assert(numCycles(T{18}{2}),1)

%!demo
%! bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! numCycles(bowtie)
%! % ans = 2
%! adj = [0 1 1; 1 0 0; 1 0 0];
%! numCycles(adj)