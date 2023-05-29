% Compute the average degree of a node in a graph, defined as
% 2 times the number of edges divided by the number of nodes
%          (every edge is counted towards the degrees twice).
%
% Inputs: adjacency matrix, nxn
% Outputs: float, the average degree, a number between 0 and max(sum(adj))
%
% Note: The average degree is related to the link density, namely:
%       link_density = ave_degree/(n-1), where n is the number of nodes
%
% Other routines used: numNodes.m, numEdges.m
% Last update, September 20, 2012

function k = averageDegree(adj)

    k = 2 * numEdges(adj) / numNodes(adj);


%!test
%!shared T
%! T = load_test_graphs();
%!assert(averageDegree(T{2}{2}),1)
%!assert(averageDegree(T{4}{2}),2+1.0/3)
%!assert(averageDegree(T{18}{2}),2)


%!demo
%! averageDegree([0 1; 1 0])
%! adj = [0 1 1; 1 0 1; 1 1 0]; % undirected 3−cycle
%! averageDegree(adj)
%! adj = [0 1 1; 1 0 0; 1 0 0];  % undirected 3−node binary tree
%! averageDegree(adj)