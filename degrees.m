% Compute the total degree, in-degree and out-degree of a graph based on the adjacency matrix;
% Note: Returns weighted degrees, if the input matrix is weighted
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: degree (1xn), in-degree (1xn) and out-degree (1xn) sequences
%
% Other routines used: isDirected.m
% Last updated, Sep 26, 2012

function [deg, indeg, outdeg] = degrees(adj)

    indeg = sum(adj);
    outdeg = sum(adj');

    if isDirected(adj)
        deg = indeg + outdeg; % total degree

    else % undirected graph: indeg=outdeg
        deg = indeg + diag(adj)'; % add self-loops twice, if any

    end


%!test
%!assert([4 4 4],degrees([0 2 1; 0 0 1; 1 1 0]))
%!shared T, deg, indeg, outdeg
%! T = load_test_graphs();
%!assert(degrees(T{1}{2}), [1 1])
%!assert(degrees(T{2}{2}), [1 1])
%!assert(degrees(T{3}{2}), [2 2])
%!assert([2 2 3 3 2 2],degrees(T{4}{2}))
%!assert([2 1 1],degrees(edgeL2adj(T{10}{2})))
%! [deg,indeg,outdeg]=degrees(edgeL2adj(T{11}{2}));
%!assert(deg,[2 1 1])
%!assert(indeg,[0 1 1])
%!assert(outdeg,[2 0 0])
%!assert(degrees(T{13}{2}), [2 2 2])
%!assert(degrees(T{14}{2}), [4 4 2])  % loops are counted twice
%!assert(degrees(T{18}{2}), [2 2 2 2])

%!demo
%! degrees([0 1 1; 1 0 1; 1 1 0])
%! [deg, indeg, outdeg] = degrees([0 1 0; 0 0 1; 1 0 0])
%! [deg, indeg, outdeg] = degrees([0 1 1; 0 0 0; 0 0 0])