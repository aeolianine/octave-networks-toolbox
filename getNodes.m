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
% GB: last updated, Jul 12 2014

function nodes = getNodes(graph,type)

if strcmp(type,'adjacency') || strcmp(type,'adjlist')
    nodes=[1:max([size(graph,1) size(graph,2)])];
 
elseif strcmp(type,'edgelist')
    nodes=unique([graph(:,1)' graph(:,2)']);
    
elseif strcmp(type,'incidence')
    nodes=[1:size(graph,1)];
else
    printf('getNodes(): "type" input can only be\n "adj" (adjacency, nxn matrix)\n, "edgelist" (mx2 or mx3 matrix)\n, "adjlist" (list, nx1 cell)\n or "inc" incidence (nxm matrix)\n')
    nodes = 'invalid graph type';
end