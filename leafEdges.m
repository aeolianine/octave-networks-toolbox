% Returns the leaf edges of the graph: edges with one adjacent edge only.
%
% Note 1: For a directed graph, leaf edges are those that "flow into" a leaf node.
% Note 2: There could be other definitions of leaves, for example: farthest away from a given root node.
% Note 3: Edges that are self-loops are not considered leaf edges.
% Note 4: Single floating disconnected edges are not considered to be leaf edges.
%
% Input: adjacency matrix, nxn
% Output: set of leaf edges: a (num edges x 2) matrix where every row contains the leaf edge nodal indices
%
% GB: last updated, Sep 23, 2012


function edges=leafEdges(adj)

adj=int8(adj>0);

lves=find(sum(adj)==1);  % same as leaf_nodes.m

edges=[];

for i=1:length(lves); edges=[edges; find(adj(:,lves(i))==1),lves(i)]; end