% Random directed graph construction
% Note 1: if p is omitted, p=0.5 is default
% Note 2: no self-loops, no double edges
%
% INPUTS:  n - number of nodes
%          p - probability, 0<=p<=1
% Output: adjacency matrix, nxn
%
% GB: last updated, Oct 21 2012

function adj = randomDirectedGraph(n,p)

adj=zeros(n); % initialize adjacency matrix

if nargin==1; p=0.5; end; % default probability


% splitting j = 1:i-1,i+1,n avoids the if statement i==j

for i=1:n

  for j=1:i-1
    if rand<=p; adj(i,j)=1; end;
  end

  
  for j=i+1:n
    if rand<=p; adj(i,j)=1; end;
  end

end

