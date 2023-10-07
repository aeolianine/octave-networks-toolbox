% Count the number of self-loops in the graph
%
% INPUT: adjacency matrix, nxn
% OUTPUT: integer, number of self-loops
%
% Note: in the adjacency matrix representation loops appear as non-zeros on the diagonal
% Last updated: Sep 20 2012

function sl = selfLoops(adj)

    sl = sum(diag(adj));


%!test
%!shared T
%! T = load_test_graphs();
%!assert(selfLoops(edgeL2adj(T{8}{2})), 1)
%!assert(selfLoops(T{14}{2}), 2)
%!assert(selfLoops(T{4}{2}),0)

%!demo
%! % no self loops
%! selfLoops([0 1; 0 0])
%! % two self loops at the same node
%! selfLoops([2 0; 0 0])
%! % three self loops
%! selfLoops([1 0 0; 0 1 0; 0 0 1])