% Checks whether a graph is simple (undirected, no self-loops, no multiple edges, no weighted edges)
%
% INPUTs: adj - adjacency matrix
% OUTPUTs: S - a Boolean variable; true (1) or false (0)
%
% Other routines used: selfLoops.m, multiEdges.m, isDirected.m
% Last updated: September 23, 2012

function S = isSimple(adj)

    S = true;

    if isDirected(adj) || selfLoops(adj) > 0 || multiEdges(adj) > 0 || length(find(adj > 1));
        S = false;
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(isSimple(T{1}{2}),false)
%!assert(isSimple(T{2}{2}),true)
%!assert(isSimple(T{3}{2}),false)
%!assert(isSimple(randomGraph(randi(5)+20,rand)),true)  % simple graph
%!assert(isSimple(edgeL2adj([1,2,2])),false)      % multi-edge
%!assert(isSimple( [1 0 0; 0 0 1; 0 1 0]),false)  % matrix with loops
%!assert(isSimple([0 1 1; 1 0 0; 0 1 0]),false)   % directed matrix

%!demo
%! % undirected binary tree 
%! isSimple([0 1 1; 1 0 0; 1 0 0])
%! % a weighted graph example
%! isSimple([0 2 1; 2 0 0; 1 0 0])
%! % directed graph example
%! isSimple([0 1 1; 0 0 0; 0 0 0])