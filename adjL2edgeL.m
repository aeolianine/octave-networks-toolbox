% Convert adjacency list to an edge list.
%
% INPUTS: adjacency list
% OUTPUTS: edge list, mx3 (m - number of edges)
%
% GB: last updated, Sep 25 2012

function el = adjL2edgeL(adjL)

el = []; % initialize edge list
for i=1:length(adjL)
    for j=1:length(adjL{i}); el=[el; i, adjL{i}(j), 1]; end
end