% Checks whether the graph is directed, using the matrix transpose function.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: boolean variable, 0 or 1
%
% Note: one-liner alternative: S=not(isSymmetric(adj));
% GB: last updated, Sep 23, 2012

function S=isDirected(adj)

S = true;
if adj==transpose(adj); S = false; end
