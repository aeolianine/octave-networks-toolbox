% Return the indices of the leaf nodes of the graph, i.e. all nodes of degree 1
% 
% Note 1: For a directed graph, leaf nodes are those with a single incoming edge
% Note 2: There could be other definitions of leaves, for example: farthest away from a given root node
% Note 3: Nodes with self-loops are not considered to be leaf nodes.
%
% Input: adjacency matrix, nxn
% Output: indices of leaf nodes
%
% GB: last updated, Sep 23, 2012

function leaves=leafNodes(adj)

adj=int8(adj>0);

leaves=find(sum(adj)==1);