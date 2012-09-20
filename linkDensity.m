##################################################################
% Computes the link density of a graph, defined as the number of edges divided by 
% number_of_nodes(number_of_nodes-1)/2 where the latter is the maximum possible number of edges.
%
% Inputs: adjacency matrix, nxn
% Outputs: link density, a float between 0 and 1
%
% Note: The graph has to be non-trivial (more than 1 node).
% Other routines used: numNodes.m, numEdges.m
% GB: last update Sep 19, 2012
##################################################################


function d=linkDensity(adj)

n = numNodes(adj);
d = 2*numEdges(adj)/(n*(n-1));