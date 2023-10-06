% Degree-preserving random rewiring.
% Note 1: Assume unweighted undirected graph.
%
% INPUTS: edge list, el (mx3) and number of rewirings, k (integer)
% OUTPUTS: rewired edge list
%
% Last updated: Sep 26, 2012

function el = rewire(el, k)

    rew = 0;

    while rew < k

        % pick two random edges
        ind = randi(length(el), 1, 2);
        edge1 = el(ind(1), :);
        edge2 = el(ind(2), :);

        if length(intersect(edge1(1:2), edge2(1:2))) > 0
            continue;
        end % the two edges cannot overlap

        % else: rewire
        if not(ismember([edge1(1), edge2(2), 1], el, 'rows')) && not(ismember([edge1(2), edge2(1), 1], el, 'rows'))

            % first possibility: (e11,e22) & (e12,e21)
            el(ind(1), :) = [edge1(1), edge2(2), 1];
            el(ind(2), :) = [edge1(2), edge2(1), 1];

            % add the symmetric equivalents
            [~, inds1] = ismember([edge1(2), edge1(1), 1], el, 'rows');
            el(inds1, :) = [edge2(2), edge1(1), 1];

            [~, inds2] = ismember([edge2(2), edge2(1), 1], el, 'rows');
            el(inds2, :) = [edge2(1), edge1(2), 1];

            rew = rew + 1;

        elseif not(ismember([edge1(1), edge2(1), 1], el, 'rows')) && not(ismember([edge1(2), edge2(2), 1], el, 'rows'))

            % second possibility: (e11,e21) & (e12,e22)
            el(ind(1), :) = [edge1(1), edge2(1), 1];
            el(ind(2), :) = [edge1(2), edge2(2), 1];

            % add the symmetric equivalents
            [~, inds1] = ismember([edge1(2), edge1(1), 1], el, 'rows');
            el(inds1, :) = [edge2(1), edge1(1), 1];

            [~, inds2] = ismember([edge2(2), edge2(1), 1], el, 'rows');
            el(inds2, :) = [edge2(2), edge1(2), 1];

            rew = rew + 1;

        else
            'rewire(): creates a double edge'
            continue
        end

    end


%!test
%! for x=1:100
%!   
%!   el = adj2edgeL(randomGraph(randi(10)+10,0.4));
%!   deg = degrees(edgeL2adj(el));
%!   rew = randi(5);
%!   eln = rewire(el,rew);
%!   degn = degrees(edgeL2adj(eln));
%!   
%!   assert(deg,degn)
%!   eq = eln(:,1:2) == el(:,1:2);
%!   preservedEdges = sum(sum(transpose(eq))==2);
%!   assert( preservedEdges >= size(eln)(1) - 4*rew )
%! 
%! end


%!demo
%! adj = randomGraph(20, 0.4);
%! elr = rewire(adj2edgeL(adj), 5);
%! adjr = edgeL2adj(elr);
%! assert(degrees(adj), degrees(adjr))
%! assert(isequal(adjr, adj), false)
