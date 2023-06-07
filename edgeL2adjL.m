% Convert an edge list to an adjacency list.
%
% INPUTS: edge list, mx3, m - number of edges
% OUTPUTS: adjacency list
%
% Note: Information about edge weights (if any) is lost.
% Last updated: September 25, 2012

function adjL = edgeL2adjL(el)

    nodes = unique([el(:, 1)' el(:, 2)']);
    adjL = cell(numel(nodes), 1);

    for e = 1:size(el, 1)
        adjL{el(e, 1)} = [adjL{el(e, 1)}, el(e, 2)];
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(edgeL2adjL(T{11}{5}),T{12}{2}')
%!assert(edgeL2adjL(sortrows(T{4}{5})),T{9}{2}')
%!assert(edgeL2adjL(T{16}{5}),T{17}{2}')


%!demo
%! edgeL2adjL([1 2 1; 1 3 1])
%! edgeL2adjL([1 2 1])