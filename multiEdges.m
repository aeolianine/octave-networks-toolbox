##################################################################
% Counts the number of multiple edges in the graph
% Multiple edges here are defined as two or more edges that have the same origin and destination nodes. 
% Note 1: This creates a natural difference in counting for undirected and directed graphs.
%
% INPUT: adjacency matrix, nxn
% OUTPUT: integer, number of multiple edges
%
% Examples: multiEdges([0 2; 2 0])=2, and multiEdges([0 0 1; 2 0 0; 0 1 0])=2
%
% Note 2: The definition of number of multi-arcs (node pairs that have multiple edges across them) 
% would be: mA = length(find(adj>1)) (normalized by 2 depending on whether the graph is directed)
%
% GB: last updated, Sep 20 2012
##################################################################


function mE=multiEdges(adj)

if adj'==adj % here use "is symmetric" as surrogate for "is directed"
  mE=sum(adj(find(adj>1)))/2;
else
  mE=sum(adj(find(adj>1)));
end