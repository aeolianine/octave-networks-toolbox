% Making an edge list (representation of a graph) symmetric,
% i.e. if [n1,n2] is in the edge list, so is [n2,n1].
%
% INPUTs: edge list, mx3
% OUTPUTs: symmetrized edge list, mx3
%
% Last updated: October 3, 2012

function el = symmetrizeEdgeL(el)

    el2 = [el(:, 1), el(:, 2)];

    for e = 1:size(el, 1)
        ind = ismember(el2, [el2(e, 2), el2(e, 1)], 'rows');
        if sum(ind) == 0; el = [el; el(e, 2), el(e, 1), el(e, 3)]; end
    end


% Alternative: Using the adjacency matrix
% adj=edgeL2adj(el);
% adj=symmetrize(adj);
% el=adj2edgeL(adj);

%!test
%! for x=1:20
%!   adj = randomDirectedGraph(randi(20)+2,rand); % create a random adjacency
%!   el = adj2edgeL(adj);
%!   if isempty(el); continue; end
%!   elsym = symmetrizeEdgeL(el);
%!   adjsym = edgeL2adj(elsym);
%!   assert(isSymmetric(adjsym),true)
%! end

%!test
%!shared T
%! T = load_test_graphs();
%!assert(sortrows(symmetrizeEdgeL(T{1}{5}))(1:2,1:2), sortrows(T{2}{5})(1:2,1:2))
%!assert(sortrows(symmetrizeEdgeL(T{6}{5}))(1:14,1:2), sortrows(T{4}{5})(1:14,1:2) )


%!demo
%! % two directed edges
%! symmetrizeEdgeL([1 2 1; 1 3 1])
%! % undirected edge: output should be the same
%! symmetrizeEdgeL([1 2 1; 2 1 1])