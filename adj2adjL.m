% Convert an adjacency graph representation to an adjacency list.
% Note 1: Valid for a general (directed, not simple) graph.
% Note 2: Edge weights (if any) get lost in the conversion.
% Note 2: Inverse function of adjL2adj.m
%
% INPUT: an adjacency matrix, nxn
% OUTPUT: cell structure for adjacency list: x{i_1}=[j_1,j_2 ...]
%
% Last updated: May 20 2023

function L = adj2adjL(adj)

    L = cell(length(adj), 1);

    for i = 1:length(adj)
        L{i} = find(adj(i, :) > 0);
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(adj2adjL(T{4}{2}),T{9}{2}')     % "bowtie" graph
%!assert(adj2adjL(T{16}{2}),T{17}{2}')   % directed 3-cycle
%%!assert(adj2adjL(edgeL2adj(T{11}{2})),T{12}{2})   % this test fails because the dimensions are flipped


%!demo
%! adj2adjL([0 1 1; 1 0 0; 1 0 0])  % undirected binary tree with 3 nodes
%! adj2adjL([0 1; 0 0])  % single directed edge