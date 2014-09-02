% Convert an edge list to an adjacency list.
% 
% INPUTS: edge list, mx3, m - number of edges
% OUTPUTS: adjacency list
%
% Note: Information about edge weights (if any) is lost.
% GB: last updated, September 25, 2012

function adjL = edgeL2adjL(el)

nodes = unique([el(:,1)' el(:,2)']);
adjL=cell(numel(nodes),1);

for e=1:size(el,1); adjL{el(e,1)}=[adjL{el(e,1)},el(e,2)]; end