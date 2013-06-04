##################################################################
% Calculate the number of independent loops (use G=m-n+c)
%      where G = num loops, m - num edges, n - num nodes, 
%                      c - number of connected components
% This is also known as the "cyclomatic number": the number of edges
% that need to be removed so that the graph doesn't have cycles.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: number of independent loops (or cyclomatic number)
%
% Other routines used: numNodes.m, numEdges.m, findConnComp.m
% GB: last updated, Oct 5 2012
##################################################################

function G = numLoops(adj)

n=numNodes(adj);
m=numEdges(adj);
comp_mat = findConnComp(adj);

G=m-n+length(comp_mat);