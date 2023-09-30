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
% Last updated: Sep 26 2014

function mE = multiEdges(adj)

    if isSymmetric(adj) % here use "is symmetric" as surrogate for "is directed"
        mE = length(find(adj > 1)) / 2;
    else
        mE = length(find(adj > 1));
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(multiEdges(T{3}{2}),1)
%!assert(multiEdges([0 2 1; 2 0 1; 1 1 0]),1)  % triangle with one double edge
%!assert(multiEdges([0 0 1; 2 0 0; 0 1 0]),1)  % directed triangle with 1 double edge
%!assert(multiEdges(randomGraph(randi(15))),0)
%!assert(multiEdges([0 0 1; 2 0 0; 0 2 0]),2)  % directed triangle with 2 double edges


%!demo
%! multiEdges([0 2; 0 0])  # one directed double edge
%! multiEdges([0 2; 2 0])  # undirected double edge
%! multiEdges([1 1 0; 1 0 0; 0 0 0])  # no multi-edges
