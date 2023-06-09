% Find the connected component to which node "i" belongs to
%
% INPUTS: adjacency matrix and index of the key node
% OUTPUTS: all node indices of the nodes in the same group
%          to which "i" belongs to (including "i")
%
% Note: Only works for undirected graphs.
% Other functions used: kneighbors.m
% Last updated: Sep 22 2012

function comp = findConnCompI(adj, i)

    neigh1 = kneighbors(adj, i, 1); % all neighbors of "i" 1 link away
    neigh1 = unique([neigh1 i]); % add i to its own component

    while 1
        len0 = length(neigh1);

        for j = 1:len0
            neigh2 = kneighbors(adj, neigh1(j), 1);
            neigh1 = unique([neigh1, neigh2]); % merge neigh1 and neigh2
        end

        if len0 == length(neigh1) % if no new neighbors found, return component
            comp = neigh1;
            return
        end

    end

end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(findConnCompI(T{5}{2},1),[1,2,3])
%!assert(findConnCompI(T{5}{2},2),[1,2,3])
%!assert(findConnCompI(T{5}{2},3),[1,2,3])
%!assert(findConnCompI(T{5}{2},4),[4,5,6])
%!assert(findConnCompI(T{5}{2},5),[4,5,6])
%!assert(findConnCompI(T{5}{2},6),[4,5,6])
%!assert(findConnCompI([0 1 0; 1 0 0; 0 0 0],1),[1,2])
%!assert(findConnCompI([0 1 0; 1 0 0; 0 0 0],2),[1,2])
%!assert(findConnCompI([0 1 0; 1 0 0; 0 0 0],3),[3])

%!demo
%! % two disconnected three−node cycles 
%! adj=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! findConnCompI(adj, 1)
%! findConnCompI(adj, 5)