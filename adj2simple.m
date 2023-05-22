% Convert an adjacency matrix of a general graph to the adjacency matrix of
%          a simple graph (symmetric, no loops, no double edges, no weights)
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: adjacency matrix (nxn) of the corresponding simple graph
%
% Other routines used: symmetrize.m
% Last updated: May 22 2023

function adj = adj2simple(adj)

    adj = adj > 0; % make all edges weight 1
    adj = symmetrize(adj);
    adj = adj - diag(diag(adj)); % clear the diagonal (selfloops)


%!test
%!assert(adj2simple(rand(6)),ones(6)-eye(6))
%!assert(adj2simple([0 2 0; 1 0 0; 1 2 0]),[0 1 1; 1 0 1; 1 1 0])
%!assert(isSymmetric(adj2simple(rand(7))),true)
