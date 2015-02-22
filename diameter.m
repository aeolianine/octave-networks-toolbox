% The longest shortest path between any two nodes nodes in the network.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: network diameter
%
% Other routines used: simpleDijkstra.m
% GB: last updated, Oct 8 2012

function diam = diameter(adj)

diam=0;
for i=1:size(adj,1)
    d=simpleDijkstra(adj,i);
    diam = max([max(d),diam]);
end