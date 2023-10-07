% This function outputs the adjacency matrix of a subgraph given
% the supergraph and the node set of the subgraph.
%
% INPUTs: adj - supergraph adjacency matrix (nxn), S - vector of subgraph node indices
% OUTPUTs: adj_sub - adjacency matrix of the subgraph (length(S) x length(S))
%
% Last update: September 23 2012

function adj_sub = subgraph(adj, S)

    adj_sub = adj(S, S);


%!test
%!shared T
%! T = load_test_graphs();
%!assert(T{13}{2},subgraph(T{4}{2},[1,2,3]))
%!assert(T{13}{2},subgraph(T{4}{2},[4,5,6]))
%!assert(T{2}{2},subgraph(T{4}{2},[4,5]))
%!assert(T{2}{2},subgraph(T{4}{2},[1,2]))
%!assert(T{2}{2},subgraph(T{4}{2},[3,4]))

%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! subgraph(bowtie, [1, 2, 3])
%! subgraph(bowtie, [3, 4, 5])