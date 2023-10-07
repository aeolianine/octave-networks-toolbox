% The signless Laplacian matrix of a graph
% Def: the sum of the diagonal degree matrix and the adjacency matrix
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: signless Laplacian matrix, nxn
%
% Last updated, Dec 6 2015

function L = signlessLaplacian(adj)

    L = diag(sum(adj)) + adj;


%!test
%!shared T
%! T = load_test_graphs();
%!assert(signlessLaplacian(T{4}{2}),[2 1 1 0 0 0; 1 2 1 0 0 0; 1 1 3 1 0 0; 0 0 1 3 1 1; 0 0 0 1 2 1; 0 0 0 1 1 2])
%!assert(signlessLaplacian(T{13}{2}),[2 1 1; 1 2 1; 1 1 2])


%!demo
%! undirected_3cycle = [0 1 1; 1 0 1; 1 1 0];
%! signlessLaplacian(undirected_3cycle)