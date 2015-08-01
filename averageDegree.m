% Compute the average degree of a node in a graph, defined as
% 2 times the number of edges divided by the number of nodes 
%          (every edge is counted towards the degrees twice).
%
% Inputs: adjacency matrix, nxn
% Outputs: float, the average degree, a number between 0 and max(sum(adj))
% 
% Note: The average degree is related to the link density, namely: 
%       link_density = ave_degree/(n-1), where n is the number of nodes
% 
% Other routines used: numNodes.m, numEdges.m
% GB: last update, September 20, 2012

function k=averageDegree(adj)

k=2*numEdges(adj)/numNodes(adj);