% Compute the link density of a graph, defined as the number of edges divided by 
% number_of_nodes(number_of_nodes-1)/2 where the latter is the maximum possible number of edges.
%
% Inputs: adjacency matrix, nxn
% Outputs: link density, a float between 0 and 1
%
% Note 1: The graph has to be non-trivial (more than 1 node).
% Note 2: Routine works for both directed and undirected graphs.
%
% Other routines used: numNodes.m, numEdges.m, isDirected.m
% GB: last update Sep 19, 2012


function d=linkDensity(adj)

n = numNodes(adj);

coeff = 2;
if isDirected(adj); coeff = 1; end

d = coeff*numEdges(adj)/(n*(n-1));