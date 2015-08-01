% Return the complement of a graph
% The complement graph has the same nodes, but edges where the original graph doesn't and vice versa.
% 
% INPUTs: adj - original graph adjacency matrix, nxn
% OUTPUTs: complement graph adjacency matrix, nxn
%
% Note 1: Assumes no multiple edges.
% Note 2: To create a complement graph without self-loops,
%         use adjC=ones(size(adj))-adj-eye(length(adj)); instead.
% GB: last updated, October 4, 2014

function adjC = graphComplement(adj)

adjC=ones(size(adj))-adj;