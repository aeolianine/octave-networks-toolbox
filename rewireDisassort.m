% Degree-preserving random rewiring.
% Every rewiring decreases the assortativity (pearson coefficient).
%
% Note 1: There are rare cases of neutral rewiring (pearson coefficient stays the same within numerical error).
% Note 2: Assume unweighted undirected graph.
%
% INPUTS: edge list, el and number of rewirings, k (integer)
% OUTPUTS: rewired edge list
%
% Other routines used: degrees.m, edgeL2adj.m
% Last updated: Sep 27 2012

function el = rewireDisassort(el, k)

    [deg, ~, ~] = degrees(edgeL2adj(el));

    rew = 0;

    while rew < k

        % pick two random edges
        ind = randi(length(el), 1, 2);
        edge1 = el(ind(1), :);
        edge2 = el(ind(2), :);

        if length(intersect(edge1(1:2), edge2(1:2))) > 0
            continue;
        end % the two edges cannot overlap

        nodes = [edge1(1) edge1(2) edge2(1) edge2(2)];
        [~, Y] = sort(deg(nodes));

        % connect nodes(Y(1))-nodes(Y(4)) and nodes(Y(2))-nodes(Y(3))
        if ismember([nodes(Y(1)), nodes(Y(4)), 1], el, 'rows') || ismember([nodes(Y(2)), nodes(Y(3)), 1], el, 'rows')
            continue;
        end

        el(ind(1), :) = [nodes(Y(1)), nodes(Y(4)), 1];
        el(ind(2), :) = [nodes(Y(2)), nodes(Y(3)), 1];

        [~, inds1] = ismember([edge1(2), edge1(1), 1], el, 'rows');
        el(inds1, :) = [nodes(Y(4)), nodes(Y(1)), 1];

        [~, inds2] = ismember([edge2(2), edge2(1), 1], el, 'rows');
        el(inds2, :) = [nodes(Y(3)), nodes(Y(2)), 1];

        rew = rew + 1;

    end


%!test
%! for x=1:100
%!   adj = [0 0; 0 0];
%!   
%!   while not(isConnected(adj)); adj = randomGraph(randi(10)+10,0.4); end
%!   el = adj2edgeL(adj);
%!   eln = rewireDisassort(el,randi(5));
%! 
%!   adjn = edgeL2adj(eln);
%!   assert(degrees(adj),degrees(adjn))
%! 
%!   assert(pearson(edgeL2adj(eln))<=(pearson(edgeL2adj(el))-10^(-7)) | (pearson(edgeL2adj(eln))-10^(-7))<=pearson(edgeL2adj(el)))
%! end


%!demo
%! adj = randomGraph(20, 0.4);
%! elr = rewireDisassort(adj2edgeL(adj), 5);
%! adjr = edgeL2adj(elr);
%! % check that the degree sequences stay the same
%! assert(degrees(adj), degrees(adjr))
%! % verify that the Pearson coefficient has decreased
%! assert(pearson(adjr)<=pearson(adj), true)