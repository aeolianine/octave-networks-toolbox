% Convert adjacency list to an edge list.
%
% INPUTS: adjacency list
% OUTPUTS: edge list, mx3 (m - number of edges)
%
% Last updated: May 20 2023

function el = adjL2edgeL(adjL)

    el = []; % initialize edge list

    for i = 1:length(adjL)

        for j = 1:length(adjL{i})
            el = [el; i, adjL{i}(j), 1];
        end

    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(adjL2edgeL(T{12}{2}),T{11}{5})           % directed 3-tree
%!assert(sortrows(adjL2edgeL(T{9}{2}))(1:14,1:2),sortrows(T{4}{5}))   % bowtie graph
%!assert(sortrows(adjL2edgeL(T{17}{2}))(1:3,1:2),T{16}{5})     % directed 3-cycle

%!demo
%! adjL2edgeL({[2 ,3], [1], [1]})
%! adjL2edgeL({[2 ,3], [], []})
%! adjL2edgeL({[2], []})