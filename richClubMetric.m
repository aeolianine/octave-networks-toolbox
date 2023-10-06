% Compute the rich club metric for a graph.
% Source: Colizza, Flammini, Serrano, Vespignani,
% "Detecting rich-club ordering in complex networks",
%                Nature Physics, vol 2, Feb 2006
%
% INPUTs: adjacency matrix, nxn, k - threshold number of links
% OUTPUTs: rich club metric
%
% Other routines used: degrees.m, subgraph.m, numEdges.m
% Last updated: October 1, 2012

function phi = richClubMetric(adj, k)

    [deg, ~, ~] = degrees(adj);

    Nk = find(deg >= k); % find the nodes with degree > k

    if isempty(Nk)
        phi = 0;
        return;
    end

    adjk = subgraph(adj, Nk);
    phi = 2 * numEdges(adjk) / (length(Nk) * (length(Nk) - 1));


%!test
%!shared T
%! T = load_test_graphs();
%!assert(richClubMetric(randomGraph(randi(5)+5,rand),12),0)
%!assert(richClubMetric(T{4}{2},2),linkDensity(T{4}{2}))
%!assert(richClubMetric(T{4}{2},3),1)
%!assert(richClubMetric(T{4}{2},4),0)

%!test
%! mat = [0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0];
%! assert(richClubMetric(mat,2),1)


%!demo
%! cycle3 = [0 1 1; 1 0 1; 1 1 0];
%! richClubMetric(cycle3, 2)
%! richClubMetric(cycle3, 1)
%! richClubMetric(cycle3, 3)
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! richClubMetric(bowtie, 1)
%! richClubMetric(bowtie, 2)
%! richClubMetric(bowtie, 3)