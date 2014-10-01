% Find the stronly connected components in a directed graph
% Source: Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms", 
%                                       SIAM Journal on Computing 1 (2): 146-160
% Wikipedia description: http://en.wikipedia.org/wiki/Tarjan's_strongly_connected_components_algorithm
% 
% Input: graph, set of nodes and edges, in adjacency list format,
%        example: L{1}=[2], L{2]=[1] is a single (1,2) edge
% Outputs: set of strongly connected components, in cell array format
%
% Other routines used: strongConnComp.m
% GB: last updated, Sep 22, 2012

function [GSCC,v] = tarjan(L)


GSCC = {};
ind = 1;                                 % node number counter 
S = [];                                  % An empty stack of nodes
for ll=1:length(L); v(ll).index = []; v(ll).lowlink = []; end  % initialize indices

for vi=1:length(L)
  if isempty(v(vi).index)
    [GSCC,S,ind,v]=strongConnComp(vi,S,ind,v,L,GSCC);  % visit new nodes only
  end
end