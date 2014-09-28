% Find the connected component to which node "i" belongs to
%
% INPUTS: adjacency matrix and index of the key node
% OUTPUTS: all node indices of the nodes in the same group 
%          to which "i" belongs to (including "i")
%
% Note: Only works for undirected graphs.
% Other functions used: kneighbors.m
% GB: last updated, Sep 22 2012


function comp=findConnCompI(adj,i)

neigh1=kneighbors(adj,i,1);   % all neighbors of "i" 1 link away
neigh1=unique([neigh1 i]);    % add i to its own component

while 1
  len0=length(neigh1);
  
  for j=1:len0
    neigh2=kneighbors(adj,neigh1(j),1);
    neigh1=unique([neigh1, neigh2]);   % merge neigh1 and neigh2
  end
  
  if len0==length(neigh1)  % if no new neighbors found, return component
    comp=neigh1;
    return
  end

end


end  % end of function