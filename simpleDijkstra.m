% Implementation of a simple version of the Dijkstra shortest path algorithm
% Returns the distances from a single vertex to all others, doesn't save the path
%
% INPUTS: adjacency matrix, adj (nxn), start node s (index between 1 and n)
% OUTPUTS: shortest path length from the start node to all other nodes, 1xn
%
% Note: Works for a weighted/directed graph.
% Last updated: September 28, 2012

function d = simpleDijkstra(adj, s)

    n = length(adj);
    d = Inf * ones(1, n); % distance s-all nodes
    d(s) = 0; % s-s distance
    T = 1:n; % node set with shortest paths not found yet

    while not(isempty(T))
        [dmin, ind] = min(d(T));

        for j = 1:length(T)

            if adj(T(ind), T(j)) > 0 && d(T(j)) > d(T(ind)) + adj(T(ind), T(j))
                d(T(j)) = d(T(ind)) + adj(T(ind), T(j));
            end

        end

        T = setdiff(T, T(ind));
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(simpleDijkstra(T{4}{2},1),[0, 1, 1, 2, 3, 3])
%!assert(simpleDijkstra(T{4}{2},3),[1, 1, 0, 1, 2, 2])

%!test
%! mat = [0 3.5 0 1; 3.5 0 1 0; 0 1 0 1.4; 1 0 1.4 0];
%! assert(simpleDijkstra(mat,1),[0, 3.4, 2.4, 1])

%!test
%!assert(simpleDijkstra(edgeL2adj(T{11}{2}),1),[0, 1, 1])
%!assert(simpleDijkstra(edgeL2adj(T{11}{2}),2),[inf, 0, inf])


%!demo
%! adj = [0 1; 0 0];
%! simpleDijkstra(adj, 1)
%! simpleDijkstra(adj, 2)
%! adj = [0 1 1; 1 0 0; 1 0 0];
%! simpleDijkstra(adj, 2)