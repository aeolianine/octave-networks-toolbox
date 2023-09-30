% Return the total number of edges given the adjacency matrix
% INPUTs: adjacency matrix, nxn
% OUTPUTs: m - total number of edges/links
%
% Note: Valid for both directed and undirected, simple or general graph
% Other routines used: selfLoops.m, isSymmetric.m
% Last updated: Sep 19, 2012

function m = numEdges(adj)

    sl = selfLoops(adj); % counting the number of self-loops

    if isSymmetric(adj) && sl == 0 % undirected simple graph
        m = sum(sum(adj)) / 2;

    elseif isSymmetric(adj) && sl > 0
        m = (sum(sum(adj)) - sl) / 2 + sl; % counting the self-loops only once

    elseif not(isSymmetric(adj)) % directed graph (not necessarily simple)
        m = sum(sum(adj));

    end


%!test
%!shared T
%! T = load_test_graphs();
%! for i=1:length(T)
%!     if strcmp(T{i}{3},'adjacency')
%!         assert( numEdges(T{i}{2}), T{i}{7} )
%!     end
%! end

%!demo
%! n = randi(100);
%! e = randi([1, n-1]);
%! adj = randomGraph(n, [], e);
%! assert(numEdges(adj), e)
%! bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! numEdges(bowtie)
