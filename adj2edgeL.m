% Convert adjacency matrix (nxn) to edge list (mx3)
%
% INPUTS: adjacency matrix: nxn
% OUTPUTS: edge list: mx3
%
% Last updated: May 20, 2023

function el = adj2edgeL(adj)

    n = length(adj); % number of nodes
    edges = find(adj > 0); % indices of all edges

    el = [];

    for e = 1:length(edges)
        [i, j] = ind2sub([n, n], edges(e)); % node indices of edge e
        el = [el; i j adj(i, j)];
    end


%!test
%!shared T
%! T = load_test_graphs();
%! for i=1:length(T)
%!     if not(strcmp( T{i}{3}, 'adjacency' )) 
%!         continue; 
%!     end
%!     edgeL1 = sortrows( adj2edgeL(T{i}{2}) );
%!     edgeL2 = sortrows( T{i}{5} );    
%!     assert(edgeL1(:,1:2), edgeL2(:,1:2))
%! end

%!demo
%! adj2edgeL([0 1 1; 1 0 0; 1 0 0])
%! adj2edgeL([0 2 1; 2 0 0; 1 0 0]) % with edge weights