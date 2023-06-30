% Check whether a graph is weighted, i.e edges have weights.
%
% INPUTS: edge list, m x 3, m: number of edges, [node 1, node 2, edge weight]
% OUTPUTS: Boolean variable, 0 or 1
%
% Last updated: Sep 23, 2012

function S = isWeighted(el)

    S = true;

    if numel(find(el(:, 3) == 1)) == size(el, 1)
        S = false;
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(isWeighted(adj2edgeL(T{2}{2})),false)
%!assert(isWeighted(adj2edgeL(T{3}{2})),true)
%!assert(isWeighted(adj2edgeL(randomGraph(randi(5)+20,rand+0.1))),false)
%!assert(isWeighted(adj2edgeL(randomDirectedGraph(randi(5)+20,rand+0.1))),false)
%!assert(isWeighted([1,2,0.5; 1,3,1.5; 1,4,1]),true)
%!assert(isWeighted([1,2,0.5; 1,3,1; 1,4,1]),true)

%!demo
%! eL =[1 2 1; 1 3 1; 2 3 1];
%! isWeighted(eL)
%! % undirected double (weighted) edge
%! isWeighted([1 2 2; 2 1 2])