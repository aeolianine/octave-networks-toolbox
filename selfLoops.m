% Count the number of self-loops in the graph
%
% INPUT: adjacency matrix, nxn
% OUTPUT: integer, number of self-loops
%
% Note: in the adjacency matrix representation loops appear as non-zeros on the diagonal
% GB: last updated, Sep 20 2012

function sl=selfLoops(adj)

sl=sum(diag(adj));