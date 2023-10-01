% Random graph construction routine.
% Note 1: Default is Erdos-Renyi graph G(n,0.5)
% Note 2: Generates undirected, simple graphs only
%
% INPUTS:  N - number of nodes
%          p - probability, 0<=p<=1
%          E - fixed number of edges; if specified, p is irrelevant
% OUTPUTS: adj - adjacency matrix of generated graph (symmetric), nxn
%
% Other routines used: numEdges.m
% Last updated: Oct 20, 2012

function adj = randomGraph(n, p, E)

    adj = zeros(n); % initialize adjacency matrix

    switch nargin % number of function arguments

        case 1% just the number of nodes, n

            % 0.5 - default probability of attachment
            for i = 1:n

                for j = i + 1:n

                    if rand <= 0.5;
                        adj(i, j) = 1;
                        adj(j, i) = 1;
                    end

                end

            end

        case 2 % the number of nodes and the probability of attachment, n, p

            for i = 1:n

                for j = i + 1:n

                    if rand <= p;
                        adj(i, j) = 1;
                        adj(j, i) = 1;
                    end

                end

            end

        case 3 % fixed number of nodes and edges, n, E

            while numEdges(adj) < E
                i = randi(n);
                j = randi(n); % pick two random nodes

                if i == j || adj(i, j) > 0;
                    continue;
                end % do not allow self-loops or double edges

                adj(i, j) = adj(i, j) + 1;
                adj(j, i) = adj(i, j);
            end

    end


%!test
%! %testing the size of the graph
%! randint = randi(20)+3;
%! assert(size(randomGraph(randint),1),randint)
%! assert(size(randomGraph(randint),2),randint)

%!test
%! % testing the default probability of attachment
%! for x=1:10
%!   randint = randi(50)+50;
%!   adj = randomGraph(randint);
%!   assert(linkDensity(adj)>0.4);
%!   assert(linkDensity(adj)<0.6);
%! end

%!test
%! % testing a random probability of attachment
%! for x=1:10
%!   p = rand;
%!   randint = randi(50)+50;
%!   adj = randomGraph(randint,p);
%!   assert(linkDensity(adj)>p-0.05);
%!   assert(linkDensity(adj)<p+0.05);
%! end

%!test
%! % testing for the number of edges, E
%! for x=1:10
%!   randint = randi(50)+50;
%!   E = randi([1,randint-1]);
%!   adj = randomGraph(randint,[],E);
%!   assert(numEdges(adj),E);
%! end


%!demo
%! adj = randomGraph(1000, 0.4);
%! numNodes(adj)
%! linkDensity(adj)