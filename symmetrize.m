% Symmetrize a non-symmetric matrix,
% i.e. returns the undirected version of a directed graph.
% Note: Where mat(i,j)~=mat(j,i), the larger (nonzero) value is chosen
%
% INPUTS: a matrix - nxn
% OUTPUT: corresponding symmetric matrix - nxn
%
% GB: last updated: October 3, 2012

function adj_sym = symmetrize(adj)

adj_sym = max(adj,transpose(adj));