% Vertex eccentricity - the maximum distance to any other vertex.
%
% Input: adjacency matrix, nxn
% Output: vector of eccentricities for all nodes, 1xn
%
% Other routines used: simpleDijkstra.m
% Last updated: Oct 10, 2012

function ec = vertexEccentricity(adj)

    n = size(adj, 1);
    ec = zeros(1, n);

    for s = 1:n;
        ec(s) = max(simpleDijkstra(adj, s));
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(vertexEccentricity(T{4}{2}),[3,3,2,2,3,3])
%!assert(vertexEccentricity(T{13}{2}),[1,1,1])
%!assert(vertexEccentricity(T{1}{2}),[1,inf])
%!assert(vertexEccentricity(edgeL2adj(T{11}{2})), [1, inf, inf])

%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! vertexEccentricity(bowtie)
%! adj = [0 1; 0 0];
%! vertexEccentricity(adj)