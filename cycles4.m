% Find cycles of length 4 in a graph;
% Note: Valid for an undirected graph only
%
% INPUTS: adjacency matrix
% OUTPUTS: number of 4-cycles, for which no edges repeat in the cycle

% Other routines used: numEdges.m, numConnTriples.m, cycles3.m
% GB: last updated January 28, 2016

function l4 = cycles4(adj)

	l4 = trace(adj^4) - 2*numEdges(adj) - 4*numConnTriples(adj) - 8*cycles3(adj);

	l4 = l4/8;