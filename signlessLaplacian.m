% The signless Laplacian matrix of a graph
% Def: the sum of the diagonal degree matrix and the adjacency matrix
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: signless Laplacian matrix, nxn
%
% Last updated, Dec 6 2015

function L = signlessLaplacian(adj)

    L = diag(sum(adj)) + adj;
