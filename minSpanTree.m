% Prim's minimal spanning tree algorithm
% Prim's alg idea:
%  start at any node, find closest neighbor and mark edges
%  for all remaining nodes, find closest to previous cluster, mark edge
%  continue until no nodes remain
%
% INPUTS: graph defined by adjacency matrix, nxn
% OUTPUTS: matrix specifying minimum spanning tree (subgraph), nxn
%
% Other routines used: isConnected.m
% GB: Oct 7, 2012

function tr = minSpanTree(adj)

% check if graph is connected:
if not(isConnected(adj)); printf('This graph is not connected. No spanning tree exists.\n'); return; end

n = length(adj); % number of nodes
tr = zeros(n);   % initialize tree

adj(find(adj==0))=inf; % set all zeros in the matrix to inf

conn_nodes = 1;        % nodes part of the min-span-tree
rem_nodes = [2:n];     % remaining nodes

while length(rem_nodes)>0
  
  [minlink]=min(min(adj(conn_nodes,rem_nodes)));
  ind=find(adj(conn_nodes,rem_nodes)==minlink);

  [ind_i,ind_j] = ind2sub([length(conn_nodes),length(rem_nodes)],ind(1));

  i=conn_nodes(ind_i); j=rem_nodes(ind_j); % gets back to adj indices
  tr(i,j)=1; tr(j,i)=1;
  conn_nodes = [conn_nodes j];
  rem_nodes = setdiff(rem_nodes,j);
  
end