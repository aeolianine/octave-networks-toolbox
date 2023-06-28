% Test whether a graph is bipartite. If so, return the two vertex sets.
% A bipartite graph is a graph for which the nodes can be split in two sets A and B,
% such that there are no edges that connect nodes within A or within B.
%
% Inputs: graph in the form of adjacency list (neighbor list, see adj2adjL.m)
% Outputs: true/false (boolean), empty set (if false) or two sets of vertices
%
% Note: This only works for undirected graphs.
% Last updated: Dec 6, 2015

function [isit, A, B] = isBipartite(L)

    isit = true; % default
    A = []; B = [];

    queue = [1]; % initialize to first vertex arbitrarily
    visited = []; % initilize to empty
    A = [1]; % put the first node on the queue in A, arbitrarily

    while not(isempty(queue))

        i = queue(1);
        visited = [visited, i];

        if length(find(A == i)) > 0

            for j = 1:length(L{i})
                B = [B, L{i}(j)];

                if length(find(visited == L{i}(j))) == 0;
                    queue = [queue, L{i}(j)];
                end

            end

        elseif length(find(B == i)) > 0

            for j = 1:length(L{i})
                A = [A, L{i}(j)];

                if length(find(visited == L{i}(j))) == 0;
                    queue = [queue, L{i}(j)];
                end

            end

        end

        queue = queue(2:length(queue)); % remove the visited node

        % if A and B overlap, return false, [],[] ....
        A = unique(A);
        B = unique(B);

        if not(isempty(intersect(A, B)));
            isit = false;
            A = [];
            B = [];
            return;
        end

    end


% Alternative for isBipartite(), without returning
% the two vertex sets.
%
% function [isit] = isBipartite(adj)
%
%    # This function uses the signlessLaplacian() function.
%    # Note that instead of checking directly for zeros,
%    # due to numeric approximation - the code checks for
%    # really small positive numbers instead.
%
%    [~,E] = eig(signlessLaplacian(adj));
%    isit = true;
%    if sum( abs(diag(E))>10^(-10) )==size(adj,1); 
%        isit = false; 
%    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(isBipartite(adj2adjL(T{4}{2})),false)
%! [isit, A, B] = isBipartite(edgeL2adjL(T{10}{2}));
%! assert(isit,true)
%! assert(A, 1)
%! assert(B, [2 3])

%! % test using the signless Laplacian
%! [~,E] = eig(signlessLaplacian(edgeL2adj(T{10}{2})));
%! isit1 = false;
%! if sum( abs(diag(E))<10^(-10) )==1
%!     isit1 = true; 
%! end
%! assert(isit, isit1)

%! even_cycle = canonicalNets(2*randi(10),'cycle');
%! [isit, A, B] = isBipartite(edgeL2adjL(even_cycle));
%! assert(isit,true)
%! assert(length(A), length(B))
%! assert(mod(A,2), ones(1,length(A)))
%! assert(mod(B,2), zeros(1,length(B)))

% test using the signless Laplacian
%! [~,E] = eig(signlessLaplacian(edgeL2adj(even_cycle)));
%! isit1 = false;
%! if sum( abs(diag(E))<10^(-10) )==1
%!     isit1 = true; 
%! end
%! assert(isit, isit1)

%! odd_cycle = canonicalNets(2*randi(10)+1,'cycle');
%! assert(isBipartite(edgeL2adjL(odd_cycle)),false)

%! % test using the signless Laplacian
%! [~,E] = eig(signlessLaplacian(edgeL2adj(odd_cycle)));
%! isit1 = false;
%! if sum( abs(diag(E))<10^(-10) )==1
%!     isit1 = true; 
%! end
%! assert(isBipartite(edgeL2adjL(odd_cycle)), isit1)

%!demo
%! % undirected binary tree with 3 nodes 
%! [is_bip, A, B] = isBipartite({[2, 3], [1], [1]})
%! % this graph contains a self−loop (2, 2)
%! [is_bip, A, B] = isBipartite({[2, 3], [1, 2], [1]})