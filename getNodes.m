% Return the list of nodes for varying graph representation types
% Inputs: graph structure (matrix or cell or struct) and type of structure (string)
%        'type' can be: 'adjacency','edgelist','adjlist' (neighbor list),'incidence' (incidence matrix)
% Note 1: only the edge list allows/returns non-consecutive node indexing
% Note 2: no build-in error check for graph structure
%
% Example representations of a directed 3-cycle: 1->2->3->1
%           'adj' - [0 1 0; 0 0 1; 1 0 0]
%           'adjlist' - {1: [2], 2: [3], 3: [1]}
%           'edgelist' - [1 2; 2 3; 3 1] or [1 2 1; 2 3 1; 3 1 1] (1 is the edge weight)
%           'inc' - [-1  0  1
%                     1 -1  0
%                     0  1 -1]
%
% Last updated: Jul 12 2014

function nodes = getNodes(graph, type)

    if strcmp(type, 'adjacency') || strcmp(type, 'adjlist')
        nodes = [1:max([size(graph, 1) size(graph, 2)])];

    elseif strcmp(type, 'edgelist')
        nodes = unique([graph(:, 1)' graph(:, 2)']);

    elseif strcmp(type, 'incidence')
        nodes = [1:size(graph, 1)];
    else
        error ('getNodes(): "type" input can only be "adjacency", "edgelist", "adjlist" or "incidence"') 
    end


%!test
%!shared T
%! T = load_test_graphs();
%! for i = 1:length(T);
%!    assert(getNodes(T{i}{2}, T{i}{3}), T{i}{4});
%! end
%! for i=1:10
%!     n = randi(100);
%!     adj = randomDirectedGraph(n);
%!     assert(getNodes(randomDirectedGraph(n),'adjacency'), 1:n)
%!     assert(getNodes(randomGraph(n),'adjacency'), 1:n)
%! end

%!test fail getNodes([], 'ergergre')