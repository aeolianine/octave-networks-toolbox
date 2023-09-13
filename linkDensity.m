% Compute the link density of a graph, defined as the number of edges divided by
% number_of_nodes(number_of_nodes-1)/2 where the latter is the maximum possible number of edges.
%
% Inputs: adjacency matrix, nxn
% Outputs: link density, a float between 0 and 1
%
% Note 1: The graph has to be non-trivial (more than 1 node).
% Note 2: Routine works for both directed and undirected graphs.
%
% Other routines used: numNodes.m, numEdges.m, isDirected.m
% Last update: Sep 19, 2012

function d = linkDensity(adj)

    n = numNodes(adj);

    coeff = 2;

    if isDirected(adj);
        coeff = 1;
    end

    d = coeff * numEdges(adj) / (n * (n - 1));


%!test
%!shared T
%! T = load_test_graphs();
%! randint = randi(101)+1;
%! assert(linkDensity(edgeL2adj(canonicalNets(randint,'tree',2))),2/randint)
%!
%! for i=1:length(T)
%!     if strcmp(T{i}{3},'adjacency')
%!         coeff = 2;
%!         if isDirected(T{i}{2}); coeff = 1; end
%!         assert( linkDensity(T{i}{2}), coeff*T{i}{7}/(T{i}{6}*(T{i}{6}-1)) )
%!     end
%! end


%!demo
%! adj = [0 1 1; 1 0 1; 1 1 0]; # undirected 3-node cycle
%! linkDensity(adj)
%! # the bowtie graph ( I>−<I ) i s a 6−node graph with 7 edges
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! linkDensity(bowtie)