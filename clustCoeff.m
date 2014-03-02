##################################################################
% Compute two clustering coefficients, based on triangle motifs count and local clustering
% C1 = number of triangle loops / number of connected triples
% C2 = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)
% Ref: M. E. J. Newman, "The structure and function of complex networks"
% Note: Valid for directed and undirected graphs
%
% INPUT: adjacency matrix, nxn
% OUTPUT: two graph average clustering coefficients (C1, C2) and clustering coefficient vector C (where mean(C) = C2)
%
% Other routines used: degrees.m, isDirected.m, kneighbors.m, numEdges.m, subgraph.m, loops3.m, numConnTriples.m
% GB, Last updated: Mar 1, 2014
% Input [in definition of C1] by Dimitris Maniadakis.
##################################################################

function [C1,C2] = clustCoeff(adj)

n = length(adj);
adj = adj>0;  % no multiple edges
[deg,~,~] = degrees(adj);
C=zeros(n,1); % initialize clustering coefficient

% multiplication change in the clust coeff formula
coeff = 2;
if isDirected(adj); coeff=1; end

for i=1:n
  
  if deg(i)==1 | deg(i)==0; C(i)=0; continue; end

  neigh=kneighbors(adj,i,1);
  edges_s=numEdges(subgraph(adj,neigh));
  
  C(i)=coeff*edges_s/(deg(i)*(deg(i)-1));

end

C1=3*loops3(adj)/(numConnTriples(adj)+2*loops3(adj));
C2=sum(C)/n;