% Check whether a graph is a tree.
% A tree is a connected graph with n nodes and (n-1) edges.
% Source: "Intro to Graph Theory" by Bela Bollobas
% 
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean variable, 0 or 1
%
% Other routines used: isConnected.m, numEdges.m, numNodes.m
% GB: last updated, Sep 24, 2012

function S=isTree(adj)

S=false;

if isConnected(adj) && numEdges(adj)==numNodes(adj)-1; S=true; end