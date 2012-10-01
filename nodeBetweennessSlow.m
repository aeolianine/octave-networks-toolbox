##################################################################
% This function returns the betweenness measure of all vertices.
% Betweenness centrality measure: number of shortest paths running through a vertex.
% Note: Valid for a general graph. Using 'number of shortest paths through a node' definition.
%
% INPUTS: adjacency or distances matrix (nxn)
% OUTPUTS: betweeness vector for all vertices (1xn)
%
% Other routines used: numNodes.m, shortestPathDP.m
% GB: Sep 28, 2012
##################################################################

function betw = nodeBetweennessSlow(adj)

n = numNodes(adj);
spaths=inf(n,n);

% calculate number of shortest paths
for k=1:n-1
  
  adjk=adj^k;

  for i=1:n
    for j=1:n
      if adjk(i,j)>0
        spaths(i,j)=min([spaths(i,j) adjk(i,j)]);
      end
    end
  end

end

betw = zeros(1,n);  

for i=1:n

    [J_st,route_st,J,route]=shortestPathDP(adj,i,i,n);

    for j=1:n
      
        if i==j; continue; end

        [J_ji,step_ind] = min(J(:,j));
        route_ji = [j, route(step_ind,j).path];
        betw(route_ji(2:length(route_ji)-1)) = betw(route_ji(2:length(route_ji)-1)) + 1/spaths(j,i);

    end

end

betw=betw/(2*nchoosek(n,2));