% Convert edge list to adjacency matrix.
%
% INPUTS: edge list: mx3, m - number of edges
% OUTPUTS: adjacency matrix nxn, n - number of nodes
%
% Note: information about nodes is lost: indices only (i1,...in) remain
% Last updated: Sep 25, 2012

function adj = edgeL2adj(el)

    nodes = sort(unique([el(:, 1) el(:, 2)])); % get all nodes, sorted
    adj = zeros(numel(nodes)); % initialize adjacency matrix

    % across all edges
    for i = 1:size(el, 1)
        adj(find(nodes == el(i, 1)), find(nodes == el(i, 2))) = el(i, 3);
    end


%!test
%!shared T
%! T = load_test_graphs();
%! for i=1:length(T)
%!     if not(strcmp( T{i}{3}, 'adjacency' )); continue; end
%!     edgeL = T{i}{5};
%!     % adding 1s to get the expected edge list dimensions right
%!     if size(edgeL)(2)==2
%!         edgeL = [edgeL ones(size(edgeL)(1),1)];
%!     end
%!     assert(T{i}{2}, edgeL2adj( edgeL ))
%! end


%!demo
%! edgeL2adj([1 2 1])
%! edgeL2adj([1 1 1; 2 3 1])