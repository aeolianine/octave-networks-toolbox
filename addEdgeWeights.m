% Add multiple edges in an edge list (sums the edge weights).
%
% INPUTS: original (non-compact) edge list
% OUTPUTS: final compact edge list (no (node1, node2) repetitions)
%
% Example: [1 2 2; 1 2 1; 4 5 1] -> [1 2 3; 4 5 1]
% Last updated: May 20 2023

function elc = addEdgeWeights(el)

    el2 = [el(:, 1), el(:, 2)]; % make the edge list searchable w/o the weights
    visited = []; % mark visited edges

    elc = [];

    for e = 1:size(el, 1)

        if sum(ismember(visited, el2(e, :), 'rows')) == 0% if not visited yet
            ind = ismember(el2, el2(e, :), 'rows');
            ind = find(ind == 1); % these are all the ocurrences of el(e,:)
            elc = [elc; el(e, 1), el(e, 2), sum(el(ind, 3))];
            visited = [visited; el2(e, :)];
        end

    end


%!test
%!assert([1 2 2; 1 3 1; 3 4 3], addEdgeWeights([1 2 1; 1 2 1; 1 3 1; 3 4 2; 3 4 1]))
%!assert([1 2 2; 2 3 4], addEdgeWeights([1 2 2; 2 3 4]))
%!assert([1 2 1; 2 1 1], addEdgeWeights([1 2 1; 2 1 1]))
%!assert([1 2 1; 2 1 1], addEdgeWeights([1 2 1; 2 1 1]))
%!assert([1 2 1; 2 1 2], addEdgeWeights([1 2 1; 2 1 1; 2 1 1]))
%!assert(addEdgeWeights([1 2 2; 1 2 1; 4 5 1]), [1 2 3; 4 5 1])

%!demo
%! addEdgeWeights([1 2 1; 1 2 0.5; 2 3 1; 3 4 1; 3 4 3])  % add the first 2 and last 2 edges
%! addEdgeWeights([1 2 1; 2 1 1; 3 4 1])   % output should be the same, no edges to add