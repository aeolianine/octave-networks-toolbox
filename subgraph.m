% This function outputs the adjacency matrix of a subgraph given 
% the supergraph and the node set of the subgraph.
%
% INPUTs: adj - supergraph adjacency matrix (nxn), S - vector of subgraph node indices
% OUTPUTs: adj_sub - adjacency matrix of the subgraph (length(S) x length(S))
%
% GB: last update, September 23 2012

function adj_sub = subgraph(adj,S)

adj_sub = adj(S,S);
