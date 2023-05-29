% Compute the average degree of neighboring nodes for every vertex.
% Note: Works for weighted degrees (graphs) also.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: average neighbor degree vector, 1xn
%
% Other routines used: degrees.m, kneighbors.m
% Last updated: Sep 28, 2012

function ave_n_deg = aveNeighborDeg(adj)

    ave_n_deg = zeros(1, length(adj)); % initialize output vector
    [deg, ~, ~] = degrees(adj);

    for i = 1:length(adj) % across all nodes

        neigh = kneighbors(adj, i, 1); % neighbors of i, one link away

        if isempty(neigh)
            ave_n_deg(i) = 0;
            continue;
        end

        ave_n_deg(i) = sum(deg(neigh)) / deg(i);

    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(aveNeighborDeg(T{13}{2}),[2 2 2])
%!assert(aveNeighborDeg(T{4}{2}),[2.5 2.5 7/3 7/3 2.5 2.5])

%!demo
%! adj = [0 1 1; 1 0 1; 1 1 0];   % undirected 3-cycle
%! aveNeighborDeg(adj)
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! aveNeighborDeg(bowtie)
%! aveNeighborDeg([0 1 1; 1 0 0; 1 0 0])