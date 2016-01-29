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
% GB, Last updated: February 6, 2015

function [C] = transitivity(adj)

C=3*cycles3(adj)/(numConnTriples(adj)+2*cycles3(adj));
