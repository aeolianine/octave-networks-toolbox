% Count the number of multiple edges in the graph.
%
% INPUT: adjacency matrix, nxn
% OUTPUT: integer, number of multiple edges
%
% Examples: multiEdges([0 2; 2 0])=1, and multiEdges([0 0 1; 2 0 0; 0 1 0])=1
%
% Note 1: The definition of number of multi-arcs/edges (node pairs
%         that have multiple edges across them) here is: 
%         mA = length(find(adj>1)); (normalized by 2 if the graph is directed).
%
% Note 2: This creates a natural difference in counting for
% undirected and directed graphs.
%
% Other routines used: isSymmetric.m
% GB: last updated, Sep 26 2014

function mE=multiEdges(adj)

if isSymmetric(adj) % here use "is symmetric" as surrogate for "is directed"
  mE=length(find(adj>1))/2;
else
  mE=length(find(adj>1));
end