% Checks whether the graph is directed, using the matrix transpose function.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: boolean variable, 0 or 1
%
% Note: one-liner alternative: S=not(isSymmetric(adj));
% Last updated: Sep 23, 2012

function S = isDirected(adj)

    S = true;

    if adj == transpose(adj)
        S = false;
    end


%!test
%!assert(isDirected(randomDirectedGraph(randi(5)+20,rand)),true)  
%!assert(isDirected(randomGraph(randi(5)+20,rand)),false)
%!shared T
%! T = load_test_graphs();
%!assert(isDirected(T{1}{2}), true)
%!assert(isDirected(T{2}{2}), false)
%!assert(isDirected(T{3}{2}), false)
%!assert(isDirected(T{16}{2}), true)

%!demo
%! isDirected([0 1 1; 1 0 0; 1 0 0])
%! isDirected([0 1 1; 0 0 0; 0 0 0])