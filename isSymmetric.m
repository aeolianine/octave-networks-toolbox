% Checks whether a matrix is symmetric (has to be square).
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: boolean variable, {0,1}
%
% GB: last update, Sep 23, 2012

function S = isSymmetric(mat)

S = false; % default
if mat == transpose(mat); S = true; end