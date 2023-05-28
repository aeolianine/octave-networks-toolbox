% Return the list of edges for varying graph representation types
% Inputs: graph structure (matrix or cell or struct) and type of structure (string)
% Outputs: edge list, mx3 matrix, where the third column is edge weight
%
% Note 1: 'type' can be: 'adjacency','edgelist','adjlist', 'incidence'
% Note 2: symmetric edges will appear twice, also in undirected graphs, (i.e. [n1,n2] and [n2,n1])
%
% Example representations of a directed triangle: 1->2->3->1
%           'adjacency' - [0 1 0; 0 0 1; 1 0 0]
%           'adjlist' - {1: [2], 2: [3], 3: [1]}
%           'edgelist' - [1 2; 2 3; 3 1] or [1 2 1; 2 3 1; 3 1 1] (1 is the edge weight)
%           'incidence' - [-1  0  1
%                           1 -1  0
%                           0  1 -1]
%
% Other routines used: adj2edgeL.m, adjL2edgeL.m, inc2edgeL.m
% Last updated: Sep 18 2012

function edges = getEdges(graph, type)

    if strcmp(type, 'adjacency')
        edges = sortrows(adj2edgeL(graph));

    elseif strcmp(type, 'edgelist')
        edges = graph; % the graph structure is the edge list

    elseif strcmp(type, 'adjlist')
        edges = sortrows(adjL2edgeL(graph));

    elseif strcmp(type, 'incidence')
        edges = sortrows(inc2edgeL(graph));
    else
        error ('getEdges(): "type" input can only be "adjacency", "edgelist", "adjlist" or "incidence"')
    end


%!test
%!shared T
%! T = load_test_graphs();
%! for i=1:length(T)
%!     edges1 = sortrows( T{i}{5} );
%!     edges2 = sortrows( getEdges(T{i}{2},T{i}{3}) ); 
%!     assert( edges1(size(edges1)(1),1:2), edges2(size(edges2)(1),1:2) )
%! end

%!test fail getEdges([], 'ergerggr')