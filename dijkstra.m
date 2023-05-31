% Dijkstra's algorithm.
%
% INPUTS: adj - adjacency matrix (nxn), s - source node, target - target node
% OUTPUTS: distance, d and path, P (from s to target)
%
% Note: if target==[], then dist and P include all distances and paths from s
% Other routines used: adj2adjL.m
% Last updated: Oct 5, 2012

function [dist, P] = dijkstra(adj, s, target)

    n = length(adj); % number of nodes
    adjL = adj2adjL(adj); % list of neighbors

    dist = inf(1, n); % initialize distances
    dist(s) = 0;

    previous = [1:n; inf(1, n)]'; % {i: inf}, i=1:n, inf -> not assigned
    S = cell(1, n); % initialize shortest path sequence

    Q = [1:n]; % all unvisited vertices, entire graph

    while length(Q) > 0 % while not empty

        % get min dist member among unvisited vertices
        [mindist, min_ind] = min(dist(Q));
        u = Q(min_ind);

        % termination condition - save "source-u" path
        S{u} = [];
        t = u;

        while not(isempty(find(previous(:, 1) == t))) % t in previous.keys():
            % insert u at the beginning of S
            S{u} = [t S{u}];
            t = previous(t, 2);
        end

        if length(target) > 0 && u == target
            dist = dist(u); P = S{u};
            return
        end

        % path book-keeping
        Q = setdiff(Q, u); % remove u from Q

        for v = 1:length(adjL{u}) % across all neighbors of u
            v = adjL{u}(v);
            alt = dist(u) + adj(u, v);

            if alt < dist(v) % update the distance and path
                dist(v) = alt;
                previous(v, 2) = u;
            end

        end

    end

    P = S;


%!test
%!shared T, d, p
%! T = load_test_graphs();
%! [d,p]=dijkstra(T{4}{2},1,5);
%!assert(d,3)
%!assert(p,[1,3,4,5])

%! [d,p]=dijkstra(T{13}{2},3,[]);
%!assert(d,[1,1,0])
%!assert(p,{[3,1],[3,2],[3]})

%! [d,p] = dijkstra(T{18}{2},3,[]);
%!assert(d,[2,1,0,1]);
%!assert(p,{[3,2,1],[3,2],[3],[3,4]})

%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! % distance and path from node 1 to node 2
%! [d,P]=dijkstra(bowtie, 1, 2)
%! [d,P]=dijkstra(bowtie, 5, 2)