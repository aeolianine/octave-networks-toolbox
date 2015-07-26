% Compute the clustering coefficient per node.
% Ci = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)
% Ref: M. E. J. Newman, "The structure and function of complex networks"
% Note: Valid for directed and undirected graphs
%
% INPUT: adjacency matrix, nxn
% OUTPUT: the average clustering coefficient (aveC) and the
%         clustering coefficient vector C per node (where mean(C) = aveC)
%
% Other routines used: degrees.m, isDirected.m, kneighbors.m,
% numEdges.m, subgraph.m
% GB, Last updated: February 7, 2015

function [aveC,C] = clustCoeff(adj)

n = length(adj);
adj = adj>0;  % no multiple edges
[deg,~,~] = degrees(adj);
C=zeros(n,1); % initialize clustering coefficient

% multiplication change in the clust coeff formula
% to reflect the edge count for directed/undirected graphs
coeff = 2;
if isDirected(adj); coeff=1; end

for i=1:n
  
  if deg(i)==1 || deg(i)==0; C(i)=0; continue; end

  neigh=kneighbors(adj,i,1);
  edges_s=numEdges(subgraph(adj,neigh));
  
  C(i)=coeff*edges_s/(deg(i)*(deg(i)-1));

end

aveC=sum(C)/n;