% Extract the giant component of a graph;
% The giant component is the largest connected component.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: giant component matrix and node indices of the giant component
%
% Other routines used: findConnComp.m, subgraph.m
% Last updated: September 22, 2012

function [GC, gc_nodes] = giantComponent(adj)

    comp = findConnComp(adj);

    L = [];

     % computing component sizes
    for k = 1:length(comp);
        L = [L, length(comp{k})];
    end

    [maxL, ind_max] = max(L);

    gc_nodes = comp{ind_max};
    GC = subgraph(adj, gc_nodes);
