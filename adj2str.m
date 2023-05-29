% Convert an adjacency matrix to a one-line string representation of a graph.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: string
%
% Note 1: The nomenclature used to construct the string is arbitrary.
%                            Here we use   .i1.j1.k1,.i2.j2.k2,....
%                            In '.i1.j1.k1,.i2.j2.k2,....',
%                            "dot" signifies new neighbor, "comma" next node
% Note 2: Edge weights are not reflected in the string representation.
% Example: [0 1 1; 0 0 0; 0 0 0] <=> .2.3,,,
%
% Other routines used: kneighbors.m
% Last updated: Sep 20 2023

function str = adj2str(adj)

    str = '';
    n = length(adj);

    for i = 1:n
        neigh = kneighbors(adj, i, 1);

        for k = 1:length(neigh)
            str = strcat(str, '.', num2str(neigh(k)));
        end

        str = strcat(str, ','); % close this node's neighbors list
    end


%!test
%!assert(adj2str(ones(3)-eye(3)),'.2.3,.1.3,.1.2,')
%!assert(adj2str(eye(3)),'.1,.2,.3,')
%!assert(adj2str([0 2; 0 0]),'.2,,')
%!shared T
%! T = load_test_graphs();
%!assert(adj2str(T{4}{2}),'.2.3,.1.3,.1.2.4,.3.5.6,.4.6,.4.5,')
%!assert(adj2str(T{16}{2}),'.2,.3,.1,')

%!demo
%! adj2str([0 1 1; 1 0 0; 1 0 0])  % undirected binary tree
%! adj2str([0 1; 0 0])  % one directed edge
%! adj2str ([1 0; 0 0])  % two disconnected nodes and a self-loop