% Return the indices of the leaf nodes of the graph, i.e. all nodes of degree 1
%
% Note 1: For a directed graph, leaf nodes are those with a single incoming edge
% Note 2: There could be other definitions of leaves, for example: farthest away from a given root node
% Note 3: Nodes with self-loops are not considered to be leaf nodes.
%
% Input: adjacency matrix, nxn
% Output: indices of leaf nodes
%
% Last updated: Sep 23, 2012

function leaves = leafNodes(adj)

    adj = int8(adj > 0);

    leaves = find(sum(adj) == 1);


%!test
%!shared T
%! T = load_test_graphs();
%!assert(leafNodes(edgeL2adj(T{10}{2})),[2,3])
%!assert(leafNodes(edgeL2adj(T{11}{2})),[2,3])
%!assert(length(leafNodes(T{13}{2})),0)
%!assert(leafNodes(T{2}{2}),[1,2])
%!assert(leafNodes(T{1}{2}),[2])
%!assert(length(leafNodes(T{4}{2})),0)
%!assert(leafNodes(edgeL2adj(T{19}{2})),[2,3,4,5])

%!demo
%! adj = [0 1 1; 1 0 0; 1 0 0];  % only ”1” is not a leaf node , because it has degree 2
%! leafNodes(adj)
%! adj = [0 1 1; 1 0 1; 1 1 0];  % a cycle has no leaf nodes
%! leafNodes(adj)
