% Find cycles of length 4 in a graph;
% Note: Valid for an undirected graph only
%
% INPUTS: adjacency matrix
% OUTPUTS: number of 4-cycles, for which no edges repeat in the cycle

% Other routines used: numEdges.m, numConnTriples.m, cycles3.m
% Last updated: January 28, 2016

function l4 = cycles4(adj)

    l4 = trace(adj^4) - 2 * numEdges(adj) - 4 * numConnTriples(adj) - 8 * cycles3(adj);

    l4 = l4 / 8;


%!test
%!shared T
%! T = load_test_graphs();
%!assert(cycles4(T{18}{2}), 1)
%!assert(cycles4(T{4}{2}),0)
%!assert(cycles4(ones(6)-eye(6)), nchoosek(6,4)*3)  % clique of size 6

%!demo
%! cycle4 = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
%! cycles4(cycle4)
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! cycles4(bowtie)