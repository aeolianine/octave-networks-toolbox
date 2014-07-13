% Convert an adjacency graph representation to an adjacency list.
% Note 1: Valid for a general (directed, not simple) graph.
% Note 2: Edge weights (if any) get lost in the conversion.
%
% INPUT: an adjacency matrix, nxn
% OUTPUT: cell structure for adjacency list: x{i_1}=[j_1,j_2 ...]
%
% GB: last updated, September 24 2012

function L = adj2adjL(adj)

L=cell(length(adj),1);

for i=1:length(adj); L{i}=find(adj(i,:)>0); end