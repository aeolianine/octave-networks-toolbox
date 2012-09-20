##################################################################
% Returns the list of edges for graph varying representation types
% Inputs: graph structure (matrix or cell or struct) and type of structure (string)
% Outputs: edge list, mx3 matrix, where the third column is edge weight
% 
% Note 1: 'type' can be: 'adj','edgelist','adjlist' (neighbor list), 'inc' (incidence matrix)
% Note 2: symmetric edges will appear twice, also in undirected graphs, (i.e. [n1,n2] and [n2,n1])
% Other routines used: adj2edgeL.m, adjL2edgeL.m, inc2edgeL.m
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


function edges = getEdges(graph,type)

if strcmp(type,'adj')
    edges=sortrows(adj2edgeL(graph));
    
elseif strcmp(type,'edgelist')
    edges=graph; % the graph structure is the edge list
    
elseif strcmp(type,'adjlist')
    edges=sortrows(adjL2edgeL(graph));
    
elseif strcmp(type,'inc')
    edges=sortrows(inc2edgeL(graph));
else
    printf('ERROR: "type" input can only be "adj" (adjacency, nxn matrix), "edgelist" (mx3 matrix)\n, "adjlist" (neighbor list, nx1 cell) and "inc" incidence (nxm matrix)\n')
end