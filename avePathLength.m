% Compute average path length for a network - the average shortest path
% Note: works for directed/undirected networks
%
% INPUTS: adjacency (or weights/distances) matrix, nxn
% OUTPUTS: average path length
%
% Other routines used: simpleDijkstra.m
% Last updated: Oct 8, 2012

function l = avePathLength(adj)

    n = size(adj, 1);

    dij = [];

    for i = 1:n
        dij = [dij; simpleDijkstra(adj, i)];
    end

    l = sum(sum(dij)) / (n^2 - n); % sum and average across everything but the diagonal


%!test
%!shared T, adj
%! T = load_test_graphs();
%!assert(avePathLength(T{4}{2}),(0+1+1+2+3+3 +0+1+2+3+3+ 0+1+2+2 +0+1+1 +0+1 +0)/15)
%!assert(avePathLength(T{13}{2}),1)
%! adj = edgeL2adj(canonicalNets(6,'line'));
%!assert(avePathLength(adj),(0+1+2+3+4+5 +0+1+2+3+4 +0+1+2+3 +0+1+2+ 0+1 +0)/15)


%!demo
%! cycle3 = [0 1 1; 1 0 1; 1 1 0]; 
%! avePathLength(cycle3)
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! avePathLength(bowtie)
