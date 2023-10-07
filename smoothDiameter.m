% A relaxed/smoothed definition of diameter: the number "d" at which
% a threshold fraction "p" of pairs of nodes are at distance at most
%                       "d". Can be non-integer using interpolation.
%
% Idea: Leskovec et al, "Graphs over Time: Densification Laws,
%                     Shrinking Diameters and Possible Explanations"
%
% Input: adjacency matrix of graph and diameter threshold, p in [0,1]
% Output: relaxed or "effective" diameter
%
% Other routines used: simpleDijkstra.m
% Last updated: Oct 8 2012

function diam = smoothDiameter(adj, p)

    n = size(adj, 1);

    dij = [];

    for i = 1:n
        dij = [dij; simpleDijkstra(adj, i)];
    end

    dij(find(dij == 0)) = inf;

    for i = 1:n - 1
        ddist(i) = length(find(dij <= i));
    end

    ddist = ddist / (n * (n - 1));

    lb = max(find(ddist <= p)); % lower bound
    ub = min(find(ddist >= p)); % upper bound

    if p == 1;
        diam = ub;
    elseif ub == lb
        diam = lb;
    elseif p < ddist(1)
        diam = 0;
    else
        % interpolation: yj = y1 + (y2-y1)*(xj-x1)/(x2-x1), where ys are diameters, xs are fractions
        diam = lb + (ub - lb) * (p - ddist(lb)) / (ddist(ub) - ddist(lb));
    end

%!test
%! adj = [0 1; 0 0];
%! while not(isConnected(adj)); adj = randomGraph(randi(10)+10,rand); end
%! assert(diameter(adj),smoothDiameter(adj,1))  % should be the same when the fraction is 1

%!test
%!shared T
%! T = load_test_graphs();
%!assert( smoothDiameter(T{13}{2},1), 1)
%!assert( smoothDiameter(T{13}{2},0.9999), 0)
%!assert( smoothDiameter(T{13}{2},0.5), 0)
%!assert( smoothDiameter(T{13}{2},0.1), 0)
%!assert( smoothDiameter(T{13}{2},0), 0)

%!assert( smoothDiameter(T{4}{2},1), 3)
%!assert( smoothDiameter(T{4}{2}, 11/15), 2)              % 7 node pairs at diameter-2
%!assert( smoothDiameter(T{4}{2}, (7/15+11/15)/2), 1.5)   % half between 1 and 2
%!assert( smoothDiameter(T{4}{2}, 7/15), 1)               % 7 node pairs at diameter-1
%!assert( smoothDiameter(T{4}{2}, 6/15), 0)
%!assert( smoothDiameter(T{4}{2}, 5/15), 0)
%!assert( smoothDiameter(T{4}{2}, 0), 0)

%!demo
%! adj = [0 1 1; 1 0 1; 1 1 0];
%! p = 1;
%! smoothDiameter(adj, p)
%! smoothDiameter(adj, 0.5)
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! smoothDiameter(bowtie, 7/15)
%! smoothDiameter(bowtie, 0.5)
%! smoothDiameter(bowtie, 9/15)
%! smoothDiameter(bowtie, 11/15)
%! smoothDiameter(bowtie, 1)
