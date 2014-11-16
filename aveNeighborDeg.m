% Compute the average degree of neighboring nodes for every vertex.
% Note: Works for weighted degrees (graphs) also.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: average neighbor degree vector, 1xn
%
% Other routines used: degrees.m, kneighbors.m
% GB: last updated, Sep 28, 2012

function ave_n_deg=aveNeighborDeg(adj)

ave_n_deg=zeros(1,length(adj));   % initialize output vector
[deg,~,~]=degrees(adj);

for i=1:length(adj)  % across all nodes
  
  neigh=kneighbors(adj,i,1);  % neighbors of i, one link away
  if isempty(neigh); ave_n_deg(i)=0; continue; end
  ave_n_deg(i)=sum(deg(neigh))/deg(i);
  
end