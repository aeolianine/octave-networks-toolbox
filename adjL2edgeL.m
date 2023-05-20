% Convert adjacency list to an edge list.
%
% INPUTS: adjacency list
% OUTPUTS: edge list, mx3 (m - number of edges)
%
% Last updated: May 20 2023

function el = adjL2edgeL(adjL)

    el = []; % initialize edge list

    for i = 1:length(adjL)

        for j = 1:length(adjL{i})
            el = [el; i, adjL{i}(j), 1];
        end

    end
