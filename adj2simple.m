% Convert an adjacency matrix of a general graph to the adjacency matrix of
%          a simple graph (symmetric, no loops, no double edges, no weights)
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: adjacency matrix (nxn) of the corresponding simple graph
%
% Other routines used: symmetrize.m
% GB: last updated, Sep 6 2014

function adj=adj2simple(adj)

adj=adj>0; % make all edges weight 1
adj = symmetrize(adj);
adj = adj - diag(diag(adj)); % clear the diagonal (selfloops)