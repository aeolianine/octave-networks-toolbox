% Return the complement of a graph
% The complement graph has the same nodes, but edges where the original graph doesn't and vice versa.
%
% INPUTs: adj - original graph adjacency matrix, nxn
% OUTPUTs: complement graph adjacency matrix, nxn
%
% Note 1: Assumes no multiple edges.
% Note 2: To create a complement graph without self-loops,
%         use adjC=ones(size(adj))-adj-eye(length(adj)); instead.
% Last updated: October 4, 2014

function adjC = graphComplement(adj)

    adjC = ones(size(adj)) - adj;


%!test
%!shared T, mat
%! T = load_test_graphs();
%! mat = [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 0 1 1; 1 1 0 1 0 0; 1 1 1 0 1 0; 1 1 1 0 0 1];
%!assert(graphComplement(T{4}{2}),mat)
%!assert(graphComplement(T{13}{2}),eye(3))
%!assert(graphComplement([0 1 1; 1 0 0; 1 0 0]), [1 0 0; 0 1 1; 0 1 1])

%!demo
%! graphComplement([0 1; 1 0])
%! graphComplement([0 1 1; 1 0 0; 1 0 0])