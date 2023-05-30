% Find cycles of length 4 in a graph and return the node indices
% Note 1: Quite basic and slow.
% Note 2: Assumes undirected graph.
%
% INPUTs: adj - adjacency matrix of graph, nxn
% OUTPUTs: 4-tuples of nodes that form 4-cycles;
%          format: {"n1-n2-n3-n4","n5-n6-n7-n8",...}
%
% Other functions used: adj2adjL.m
% Last updated: Oct 5 2012

function l4 = cycle4nodes(adj)

    n = size(adj, 1); % number of nodes
    L = adj2adjL(adj); % adjacency list or list of neighbors

    l4 = {}; % initialize loops of size 4

    for i = 1:n - 1

        for j = i + 1:n

            int = intersect(L{i}, L{j});
            int = setdiff(int, [i j]);

            if length(int) >= 2
                % enumerate pairs in the intersection
                for ii = 1:length(int) - 1

                    for jj = ii + 1:length(int)
                        loop4 = sort([i, j, int(ii), int(jj)]);
                        loop4 = strcat(num2str(loop4(1)), '-', num2str(loop4(2)), '-', num2str(loop4(3)), '-', num2str(loop4(4)));

                        if sum(ismember(l4, loop4)) > 0;
                            continue;
                        end

                        l4{length(l4) + 1} = loop4;

                    end

                end

            end

        end

    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(cycle4nodes(T{18}{2}), {'1-2-3-4'})
%!assert(cycle4nodes(T{4}{2}),{})
%!assert(cycle4nodes(ones(4)-eye(4)),{'1-2-3-4'})  % clique of size 4
%!assert(length(cycle4nodes(ones(6)-eye(6))),nchoosek(6,4))  % clique of size 6

%!demo
%! cycle_four_nodes = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]; 
%! cycle4nodes(cycle_four_nodes)