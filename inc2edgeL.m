% Convert an incidence matrix to an edge list.
% 
% Inputs: inc - incidence matrix nxm (number of nodes x number of edges)
% Outputs: edge list - mx3, m x (node 1, node 2, edge weight)
%
% Example: [-1; 1] <=> [1,2,1], one directed (1->2) edge
% GB: last updated, Sep 25 2012

function el = inc2edgeL(inc)

m = size(inc,2); % number of edges
el = zeros(m,3); % initialize edge list [n1, n2, weight]

for e=1:m
    ind_m1 = find(inc(:,e)==-1);
    ind_p1 = find(inc(:,e)==1);
    
    if numel(ind_m1)==0 && numel(ind_p1)==1  % undirected, self-loop
        el(e,:) = [ind_p1 ind_p1 1];  
        
    elseif numel(ind_m1)==0 && numel(ind_p1)==2 % undirected
        el(e,:) = [ind_p1(1) ind_p1(2) 1];
        el=[el; ind_p1(2) ind_p1(1) 1];
        
    elseif numel(ind_m1)==1 && numel(ind_p1)==1 % directed
        el(e,:) = [ind_m1 ind_p1 1];
        
    end
end