% Check whether a graph is complete, i.e. whether every node is linked to every other node.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean variable, true/false
%
% Note: Only defined for unweighted graphs.
% GB: last updated, Sep 23, 2012

function S=isComplete(adj)

S=false; % default

adj=adj>0;  % remove weights
n=length(adj);

% all degrees "n-1" or "n" or w/ n selfloops
if sum(adj)==ones(1,n)*(n-1) || sum(adj)==ones(1,n)*n; S=true; end