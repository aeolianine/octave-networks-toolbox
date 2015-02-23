% Checks whether a graph is simple (undirected, no self-loops, no multiple edges, no weighted edges)
%
% INPUTs: adj - adjacency matrix
% OUTPUTs: S - a Boolean variable; true (1) or false (0)
%
% Other routines used: selfLoops.m, multiEdges.m, isDirected.m
% GB: last updated, September 23, 2012

function S = isSimple(adj)

S=true;

if isDirected(adj) || selfLoops(adj)>0 || multiEdges(adj)>0 || length(find(adj>1)); S=false; end