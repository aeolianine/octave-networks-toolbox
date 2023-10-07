% Symmetrize a non-symmetric matrix,
% i.e. returns the undirected version of a directed graph.
% Note: Where mat(i,j)~=mat(j,i), the larger (nonzero) value is chosen
%
% INPUTS: a matrix - nxn
% OUTPUT: corresponding symmetric matrix - nxn
%
% Last updated: October 3, 2012

function adj_sym = symmetrize(adj)

    adj_sym = max(adj, transpose(adj));


%!test
%! for i=1:20
%!  adj = randomDirectedGraph(randi(10)+3,rand);
%!  assert(isSymmetric(symmetrize(adj)),true)
%! end

%!test
%!shared T
%! T = load_test_graphs();
%!assert(symmetrize(T{1}{2}),T{2}{2})
%!assert(symmetrize(edgeL2adj(T{11}{2})),edgeL2adj(T{10}{2}))
%!assert(symmetrize(T{16}{2}),T{13}{2})


%!demo
%! symmetrize([0 1; 0 0])