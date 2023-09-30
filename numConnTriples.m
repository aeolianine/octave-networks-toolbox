% Count the number of connected triples in a graph.
% Note: works for undirected graphs only
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: integer - number of connected triples
%
% Other routines used: kneighbors.m, cycles3.m
% Last updated: October 4, 2012

function c = numConnTriples(adj)

    c = 0; % initialize

    for i = 1:length(adj)
        neigh = kneighbors(adj, i, 1);

        if length(neigh) < 2
            continue;
        end % handle leaves, no triple here

        c = c + nchoosek(length(neigh), 2);
    end

    c = c - 2 * cycles3(adj); % due to the symmetry triangles repeat 3 times
    %                        in the nchoosek count

    % alternative
    % def numConnTriples(L):
    %
    % % input: adjacency list
    % % outputs: number of connected triples
    %
    % c=0;      % initialize number of connected triples
    %
    % for i=1:length(L)
    %   neigh = L{i}
    %   if length(neigh)<2: continue; end
    %   c = c + nchoosek(length(neigh),2);
    % end
    %
    % c = c - 2*cycles3(adjL2adj(L));


%!test
%!shared T
%! T = load_test_graphs();
%!assert(numConnTriples(T{4}{2}),6)
%!assert(numConnTriples(T{13}{2}),1)
%!assert(numConnTriples(edgeL2adj(T{10}{2})),1)
%!assert(numConnTriples(T{2}{2}),0)
%!assert(numConnTriples(T{18}{2}),4)

%!demo
%! cycle3 = [0 1 1; 1 0 1; 1 1 0];
%! numConnTriples(cycle3)
%! % ans = 1
%! adj = [0 1 0; 1 0 0; 0 0 0];
%! numConnTriples(adj)
%! % ans = 0
%! bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! numConnTriples(bowtie)


