% Sort nodes by degree, and when there's equality, by maximum neighbor degree
% Ideas from Guo, Chen, Zhou, "Fingerprint for Network Topologies"
%
% INPUTS: adjacency matrix, 0s and 1s, nxn
% OUTPUTS: sorted (decreasing) sequence of nodal indices (nx1)
%
% Note: Works for undirected graphs only.
% Other routines used: degrees.m, kneighbors.m
% Last updated, Nov 24, 2014

function I = sortNodesByMaxNeighborDegree(adj)

    [deg, ~, ~] = degrees(adj); % compute all degrees, use "deg" assuming symmetry

    degmat = zeros(size(adj, 1), 2); % a nx2 matrix

    for x = 1:size(adj, 1)% across all nodes
        degmat(x, 1) = deg(x);

        if deg(x) == 0
            degmat(x, 2) = 0;
            continue
        end

        nei_inds = kneighbors(adj, x, 1);
        degmat(x, 2) = max(deg(nei_inds));

    end

    [sortmat, I] = sortrows(degmat); % increasing sorting
    I = I(end:-1:1); % decreasing sorting
