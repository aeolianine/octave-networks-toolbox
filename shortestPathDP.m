% Shortest path algorithm using dynamic programming.
% Note 1: Valid for directed/undirected network.
% Note 2: if links have weights, they are treated as distances.
% Source: D. P. Bertsekas, Dynamic Programming and Optimal Control, 
%                             Athena Scientific, 2005 (3rd edition)
%
% INPUTs: L - (cost/path lengths matrix), s - (start/source node), 
%                                       t - (end/destination node)
%                                 steps - number of arcs allowable
% OUTPUTS: 
%       route - sequence of nodes on optimal path, at current stage
%       route(k,i).path - best route from "i" to destination "t" in "k" steps
%       route_st - best route from "s" to "t"
%       J_st - optimal cost function (path length) from "s" to "t"
%       J(1,i) - distance from node "i" to "t" in "k" steps
%
% GB: last updated, Oct 5 2012

function [J_st,route_st,J,route]=shortestPathDP(L,s,t,steps)

n = size(L,2);

L(find(L==0))=Inf;  % make all zero distances equal to infinity

for i=1:n
  J(steps,i) = L(i,t); 
  route(steps,i).path = [t];
end

% find min for every i: Jk(i)=min_j(L(i,j)+Jk+1(j))
for p=1:steps-1
  k=steps-p; % recurse backwards
  
  for i=1:n
    [J(k,i),ind_j] =  min(L(i,:)+J(k+1,:));
    route(k,i).path = [ind_j, route(k+1,ind_j).path];
  end
  
end



[J_st,step_ind] = min(J(:,s));               % the shortest path (min cost) from s to t
route_st = [s, route(step_ind,s).path];      % the shortest path route from s to t
J=J(sort(1:min([n,steps]),'descend'),:);
route=route(sort(1:min([n,steps]),'descend'),:);