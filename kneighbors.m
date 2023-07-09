% Finds the number of k-neighbors (k links away) for every node
%
% INPUTS: adjacency matrix (nxn), start node index, k - number of links
% OUTPUTS: vector of k-neighbors indices
%
% Last updated: Oct 7 2012

function kneigh = kneighbors(adj, ind, k)

    adjk = adj;

    for i = 1:k - 1
        adjk = adjk * adj;
    end;

    kneigh = find(adjk(ind, :) > 0);


%!test
%!shared T
%! T = load_test_graphs();
%!assert(kneighbors(T{4}{2},1,3),[1 2 3 4 5 6])
%!assert(kneighbors(T{4}{2},3,1),[1 2 4])
%!assert(kneighbors(T{13}{2},2,1),[1,3])
%!assert(kneighbors(T{13}{2},1,2),[1,2,3])

%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! kneighbors(bowtie, 1, 1)
%! kneighbors(bowtie, 1, 2)
%! kneighbors(bowtie, 1, 3)
%! kneighbors(bowtie, 4, 2)

