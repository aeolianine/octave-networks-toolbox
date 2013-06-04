##################################################################
% Returns the list of nodes for varying graph representation types
% Inputs: graph structure (matrix or cell or struct) and type of structure (string)
%        'type' can be: 'adj','edgelist','adjlist' (neighbor list),'inc' (incidence matrix)
% Note 1: only the edge list allows/returns non-consecutive node indexing
% Note 2: no build-in error check for graph structure
%
% Example representations of a directed triangle: 1->2->3->1
%           'adj' - [0 1 0; 0 0 1; 1 0 0]
%           'adjlist' - {1: [2], 2: [3], 3: [1]}
%           'edgelist' - [1 2; 2 3; 3 1] or [1 2 1; 2 3 1; 3 1 1] (1 is the edge weight)
%           'inc' - [-1  0  1
%                     1 -1  0
%                     0  1 -1]
%
% GB: last updated, Sep 18 2012
##################################################################

function nodes = getNodes(graph,type)

if strcmp(type,'adj') | strcmp(type,'adjlist')
    nodes=[1:max([size(graph,1) size(graph,2)])];
 
elseif strcmp(type,'edgelist')
    nodes=unique([graph(:,1)' graph(:,2)']);
    
elseif strcmp(type,'inc')
    nodes=[1:size(graph,1)];
else
    printf('ERROR: "type" input can only be "adj" (adjacency, nxn matrix), "edgelist" (mx2 or mx3 matrix)\n, "adjlist" (neighbor list, nx1 cell) and "inc" incidence (nxm matrix)\n')
    
end