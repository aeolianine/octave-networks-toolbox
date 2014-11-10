% Compute the total degree, in-degree and out-degree of a graph based on the adjacency matrix; 
% Note: Returns weighted degrees, if the input matrix is weighted
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: degree (1xn), in-degree (1xn) and out-degree (1xn) sequences
%
% Other routines used: isDirected.m
% GB: last updated, Sep 26, 2012

function [deg,indeg,outdeg]=degrees(adj)

indeg = sum(adj);
outdeg = sum(adj');

if isDirected(adj)
  deg = indeg + outdeg; % total degree

else   % undirected graph: indeg=outdeg
  deg = indeg + diag(adj)';  % add self-loops twice, if any

end