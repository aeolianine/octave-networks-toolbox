##################################################################
% Build a random modular graph, given number of modules, and link densities
%
% INPUTs: number of nodes (n), number of modules (c), total link density (p),
%         and ratio of nodal degree to nodes within the same module
%         to the degree to nodes in other modules (r)
% OUTPUTs: adjacency matrix, modules to which the nodes are assigned
%
% GB: last updated, November 4, 2012
##################################################################

function [adj, modules] = randomModularGraph(n,c,p,r)

% n - number of nodes
% c - number of clusters/modules
% p - overall probability of attachment
% r - proportion of links within modules


% assign nodes to modules: 1 -> n/c, n/c+1 -> 2n/c, ... , (c-1)n/c -> c(n/c);
modules={};
for k=1:c; modules{k}=round((k-1)*n/c+1):round(k*n/c); end

adj=zeros(n); % initialize adjacency matrix

% DERIVATION of probabilities
% k_in/k_out = r, k_in + k_out = k = p(n-1)
% => (1/r)k_in + k_in = p(n-1), => k_in = rp(n-1)/(r+1), k_out = p(n-1)/(r+1)
% k_in = p_in*(n/c-1) => p_in = rpc(n-1)/((r+1)(n-c))
% k_out = p_out*(n-n/c) => p_out = pc(n-1)/(n(r+1)(c-1))


p_in = r*p*c*(n-1)/((r+1)*(n-c));
p_out = p*c*(n-1)/(n*(r+1)*(c-1));


for i=1:n
  for j=i+1:n
        
      module_i=ceil(c*i/n);   % the module to which i belongs to
      module_j=ceil(c*j/n);   % the module to which j belongs to
        
      if module_i==module_j
        
        % prob of attachment within module
        if rand<p_in; adj(i,j)=1; adj(j,i)=1; end
        
      else
      
        % prob of attachment across modules
        if rand<p_out; adj(i,j)=1; adj(j,i)=1; end
        
      end
  end
end