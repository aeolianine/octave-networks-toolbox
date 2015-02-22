% Compute average path length for a network - the average shortest path
% Note: works for directed/undirected networks 
%
% INPUTS: adjacency (or weights/distances) matrix, nxn
% OUTPUTS: average path length
%
% Other routines used: simpleDijkstra.m 
% GB: Oct 8, 2012

function l = avePathLength(adj)

n=size(adj,1);

dij = [];

for i=1:n; dij=[dij; simpleDijkstra(adj,i) ]; end

l = sum(sum(dij))/(n^2-n); % sum and average across everything but the diagonal