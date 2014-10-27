% Check whether a graph is weighted, i.e edges have weights.
% 
% INPUTS: edge list, m x 3, m: number of edges, [node 1, node 2, edge weight]
% OUTPUTS: Boolean variable, 0 or 1
%
% GB: last updated, Sep 23, 2012

function S=isWeighted(el)

S=true;

if numel( find(el(:,3)==1) ) == size(el,1); S=false; end