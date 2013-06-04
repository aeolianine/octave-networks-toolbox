##################################################################
% Returns the complement of a graph
% The complement graph has the same nodes, but edges where the original graph doesn't and vice versa.
% 
% INPUTs: adj - original graph adjacency matrix, nxn
% OUTPUTs: complement graph adjacency matrix, nxn
%
% Note: Assumes no multiple edges
% GB: last updated, September 23, 2012
##################################################################

function adjC = graphComplement(adj)

adjC=ones(size(adj))-adj;