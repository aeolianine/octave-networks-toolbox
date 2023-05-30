% Compute the closeness centrality for every vertex: 1/sum(dist to all other nodes)
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: vector of closeness centralities, nx1
%
% Source: social networks literature (example: Wasserman, Faust, "Social Networks Analysis")
% Other routines used: simpleDijkstra.m
% Last updated: Sep 28, 2012

function C = closeness(adj)

    C = zeros(length(adj), 1); % initialize closeness vector

    for i = 1:length(adj)
        C(i) = 1 / sum(simpleDijkstra(adj, i));
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(closeness(T{4}{2})',[1/(1+1+2+3+3), 1/(1+1+2+3+3), 1/(1+1+1+2+2), 1/(1+1+1+2+2), 1/(1+1+2+3+3), 1/(1+1+2+3+3)])
%!assert(closeness([0 1 1; 1 0 0; 1 0 0]),[0.5 1/3 1/3]')
%!assert(closeness(T{13}{2}),[1/(1+1), 1/(1+1), 1/(1+1)]')

%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! closeness(bowtie)
%! adj = [0 1 1; 1 0 1; 1 1 0]; 
%! closeness(adj)