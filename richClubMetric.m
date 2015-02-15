% Compute the rich club metric for a graph.
% Source: Colizza, Flammini, Serrano, Vespignani, 
% "Detecting rich-club ordering in complex networks", 
%                Nature Physics, vol 2, Feb 2006
% 
% INPUTs: adjacency matrix, nxn, k - threshold number of links
% OUTPUTs: rich club metric
%
% Other routines used: degrees.m, subgraph.m, numEdges.m
% GB: last updated, October 1, 2012

function phi=richClubMetric(adj,k)

[deg,~,~]=degrees(adj);

Nk=find(deg>=k);       % find the nodes with degree > k
if isempty(Nk); phi = 0; return; end

adjk=subgraph(adj,Nk);
phi=2*numEdges(adjk)/(length(Nk)*(length(Nk)-1));