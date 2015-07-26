% Implements the network growth model from: Fabrikant et al, 
%                  "Heuristically Optimized Trade-offs: A New Paradigm
%                                       for Power Laws in the Internet"
% Note: Assumes the central point (root) to be the one closest to (0,0)
%
% INPUTS: n - number of points, parameter alpha, [0,inf), 
%                     plt='on'/'off', if [], then 'off' is default
% OUTPUTS: generated network (adjacency matrix) and plot [optional]
%
% GB: last updated: November 14, 2012


function [adj,p]=fabrikantModel(n,alpha,plt)


% create random point (polar) coordinates
p = zeros(n,2);
rad = rand(1,n);
theta = 2*pi*rand(1,n);
for i=1:n; p(i,:) = rad(i)*[cos(theta(i)),sin(theta(i))]; end
norms = [];
for pp=1:length(p); norms = [norms; norm(p(pp,:))]; end
[~,Y]=sort(norms);
p = p(Y,:);
clear norms rad theta

h=zeros(n,1); % initialize centrality function
adj=zeros(n); % initialize adjacency matrix

adjL = pdist2(p,p,'euclidean');  % compute all point-to-point distances

% compute centrality function based from p(1,:)
h = adjL(1,:);


for i=2:n   % a new point arrives at each iteration

  % compute min weight across all existing points
  d=[]; 
  for j=1:i-1; d=[d; alpha*norm(p(i,:)-p(j,:)) + h(j)]; end
  [~,indmin]=min(d);
  adj(i,indmin)=1; adj(indmin,i)=1;
  
  
end

if nargin<=2 || not(strcmp(plt,'on')); return; end

set(gcf,'color',[1,1,1])
for i=1:n
  axis off
  %drawnow
  %plot(p(i,1),p(i,2),'.','Color',[0.5,0.5,0.5]); hold off; hold on;  % plot the new point
  for j=i+1:n
    if adj(i,j)>0
      line([p(i,1) p(j,1)],[p(i,2) p(j,2)],'Color',[0.5,0.5,0.5])
      hold off; hold on;
    end
  end
end





function L = pdist2(x,y,distance);

% L is the matrix of distances between 
% all members of x to all members of y
% distance type can be defined differently,
% but for the purposes of this routine it
% is 'euclidean'

if not(strcmp(distance,'euclidean'));
  print('no other distance than euclidean implemented yet\n');
  return
end

L = zeros(size(x,1),size(y,1));

for xx=1:size(x,1)
  for yy=1:size(y,1)
    L(xx,yy) = norm(x(xx,:)-y(yy,:));
  end
end