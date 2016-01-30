% Calculate the number of independent cycles (use L = m-n+c)
%      where L = num cycles, m - num edges, n - num nodes, 
%                      c - number of connected components
% This is also known as the "cyclomatic number": the number of edges
% that needs to be removed so that the graph doesn't have cycles.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: number of independent cycles (or cyclomatic number)
%
% Other routines used: numNodes.m, numEdges.m, findConnComp.m
% GB: last updated, Oct 5 2012

function L = numCycles(adj)

n=numNodes(adj);
m=numEdges(adj);
comp_mat = findConnComp(adj);

L=m-n+length(comp_mat);