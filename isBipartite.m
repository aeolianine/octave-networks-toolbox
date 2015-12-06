% Test whether a graph is bipartite. If so, return the two vertex sets.
% A bipartite graph is a graph for which the nodes can be split in two sets A and B, 
% such that there are no edges that connect nodes within A or within B.
% 
% Inputs: graph in the form of adjacency list (neighbor list, see adj2adjL.m)
% Outputs: true/false (boolean), empty set (if false) or two sets of vertices
%
% Note: This only works for undirected graphs.
% GB: last updated, Dec 6, 2015

function [isit,A,B]=isBipartite(L)

isit=true;   % default
A=[]; B=[];

queue=[1];   % initialize to first vertex arbitrarily
visited=[];  % initilize to empty
A=[1];       % put the first node on the queue in A, arbitrarily

while not(isempty(queue))
  
  i=queue(1);
  visited=[visited, i];
  
  if length(find(A==i))>0
    for j=1:length(L{i})
      B=[B,L{i}(j)];
      if length(find(visited==L{i}(j)))==0; queue=[queue, L{i}(j)]; end
    end
      
  elseif length(find(B==i))>0
   
    for j=1:length(L{i})
      A=[A,L{i}(j)];
      if length(find(visited==L{i}(j)))==0; queue=[queue, L{i}(j)]; end
    end
  
  end
  
  queue=queue(2:length(queue)); % remove the visited node
  
  % if A and B overlap, return false, [],[] ....
  A=unique(A); B=unique(B);
  if not(isempty(intersect(A,B))); isit=false; A=[]; B=[]; return; end
  % ............................................
  
end


% Alternative for isBipartite(), without returning 
% the two vertex sets.
%
% function [isit] = isBipartite(adj)
%
%    # This function uses the signlessLaplacian() function.
%    # Note that instead of checking directly for zeros,
%    # due to numeric approximation - the code checks for
%    # really small positive numbers instead.
%
%    [~,E] = eig(signlessLaplacian(adj));
%    isit = true;
%    if sum( abs(diag(E))>10^(-10) )==size(adj,1); isit = false; end

