% "Master equation" growth model, as presented in 
%                   "Evolution of Networks" by Dorogovtsev, Mendez
% Note: probability of attachment: (q(i)+ma)/((1+a)mt), 
%             q(i)-indegree of i, a=const, t - time step (# nodes)
%
% INPUTS: number of nodes n, m - # links to add at each step, a=constant
% OUTPUTS: adjacency matrix, nxn
%
% Other routines used: weightedRandomSample.m
% GB: last updated, Nov 11, 2012

function adj=masterEquationGrowthModel(n,m,a)

adj=zeros(n); adj(1,2)=2; adj(2,1)=2; % initial condition
vertices = 2; 
if nargin==2 || a==[]; a = 2; end      % pick a constant

while vertices < n
  
  t = vertices;
  q=sum(adj); % in-degrees
    
  % compute the probability of attachment
  pk = zeros(1,t);
  for k=1:t; pk(k)=(q(k)+m*a)/((1+a)*m*t); end

  r = weightedRandomSample(m,[1:t],pk/sum(pk));
  if length(unique(r))~=length(r)
    r = weightedRandomSample(m,[1:t],pk/sum(pk));
  end
  
  vertices=vertices+1;   % add vertex

  % add m links
  for node=1:length(r)
    adj(vertices,r(node))=1;
    adj(r(node),vertices)=1;
  end
    
  adj(vertices,vertices)=1;  % assure non-zero probability of attachment

end

adj = adj>0;              % no multi-edges
adj=adj-diag(diag(adj));  % remove self-loops