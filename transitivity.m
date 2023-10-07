% Calculate the transitivity.
% C = number of 3-cycles / number of connected triples
% Ref: M. E. J. Newman, "The structure and function of complex networks"
% Note: Valid for directed and undirected graphs
%
% INPUT: adjacency matrix, nxn
% OUTPUT: The transitivity, C
%
% Other routines used: cycles3.m, numConnTriples.m
% Input/corrections by Dimitris Maniadakis.
% Last updated: February 6, 2015

function [C] = transitivity(adj)

    C = 3 * cycles3(adj) / (numConnTriples(adj) + 2 * cycles3(adj));


%!test
%!shared T
%! T = load_test_graphs();
%!assert(clustCoeff(T{13}{2}),1)
%!assert(clustCoeff(edgeL2adj(T{10}{2})),0)
%!assert(clustCoeff(edgeL2adj(canonicalNets(randi(10)+5,'tree',2))),0)
%!assert(transitivity(T{4}{2}),0.6)


%!demo
%! adj = [0 1 1; 1 0 1; 1 1 0];
%! transitivity(adj)
%! adj = [0 1 1; 1 0 0; 1 0 0];
%! transitivity(adj)
%! adj = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! transitivity(adj)