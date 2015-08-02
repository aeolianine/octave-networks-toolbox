% Compute the closeness centrality for every vertex: 1/sum(dist to all other nodes)
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: vector of closeness centralities, nx1
%
% Source: social networks literature (example: Wasserman, Faust, "Social Networks Analysis")
% Other routines used: simpleDijkstra.m 
% GB: last updated, Sep 28, 2012

function C=closeness(adj)

C=zeros(length(adj),1);  % initialize closeness vector

for i=1:length(adj); C(i)=1/sum( simpleDijkstra(adj,i) ); end