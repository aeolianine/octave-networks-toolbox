% The ith component of the eigenvector corresponding to the greatest
% eigenvalue gives the centrality score of the ith node in the network.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: eigen(-centrality) vector, nx1
%
% Last updated: Sep 29, 2012

function x = eigenCentrality(adj)

    [V, D] = eig(adj);
    [~, ind] = max(diag(D));
    x = V(:, ind);


%!test
%!shared T, v, ec, adj
%! T = load_test_graphs();
%! [v,~]=eig([0 1 1; 1 0 1; 1 1 0]);
%!assert(eigenCentrality([0 1 1; 1 0 1; 1 1 0]),v(:,3))   % "3" is the number of nodes

%! [v,~]=eig(T{4}{2});
%!assert(eigenCentrality( T{4}{2} ),v(:,size(T{4}{2},1)))

%! [v,~]=eig(T{13}{2});
%!assert(eigenCentrality( T{13}{2} ),v(:,size(T{13}{2},1)))

%! [v,~]=eig(T{18}{2});
%! ec = v(:,size(T{18}{2},1));
%!assert(eigenCentrality( T{18}{2} ), ec)
%!assert(norm( ec(1)*ones(length(ec),1) - ec) < 1*e^(-20))

%! adj = edgeL2adj(canonicalNets(randi(10)+2, 'cycle'));
%! [v,~]=eig(adj);
%! ec = v(:,size(adj,1));
%!assert(eigenCentrality( adj ), ec)
%!assert(norm( ec(1)*ones(length(ec),1) - ec) < 1*e^(-20))


%!demo
%! eigenCentrality([0 1 1; 1 0 1; 1 1 0])
%! adj = [0 1 1; 1 0 0; 1 0 0]; 
%! eigenCentrality(adj)