% Vertex eccentricity - the maximum distance to any other vertex.
%
% Input: adjacency matrix, nxn
% Output: vector of eccentricities for all nodes, 1xn
%
% Other routines used: simpleDijkstra.m
% GB: last updated, Oct 10, 2012

function ec=vertexEccentricity(adj)

n=size(adj,1);
ec=zeros(1,n);

for s=1:n; ec(s)=max( simpleDijkstra(adj,s) ); end