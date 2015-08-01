% Check whether a graph is regular, i.e. whether every node has the same degree.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean, 0 or 1
%
% Note: Defined for unweighted graphs only.
% GB: last updated, Sep 23, 2012

function S=isRegular(adj)

S=false;

degs=sum(adj>0); % remove weights and sum columns

if degs == degs(1)*ones(size(degs)); S = true; end