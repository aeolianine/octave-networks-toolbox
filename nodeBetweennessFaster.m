##################################################################
% Betweenness centrality measure: number of shortest paths running through a vertex.
% Compute for all vertices, using Dijkstra's algorithm, and the number
%                        of shortest paths through a node' definition.
% Note 1: The graph has to be connected.
% Note 2: Bug: currently the routine does not return the correct
%       betweenness values for all nodes if the graph contains an even
%       cycle. That is because an even cycle results in multiple shortest
%       paths between some of the nodes, and this function chooses one
%       shortest path for every pair of nodes. Fix TBA.
%
% INPUTS: adjacency or distances matrix, nxn
% OUTPUTS: betweeness vector for all vertices (1xn)
%
% Other routines used: dijkstra.m
% GB: June 4 2013
##################################################################

function betw = nodeBetweennessFaster(adj)

n = length(adj);
spaths=inf(n,n);
adjk = adj;

% calculate number of shortest paths
for k=1:n-1
  
  for i=1:n
    for j=1:n
      if adjk(i,j)>0; spaths(i,j)=min([spaths(i,j) adjk(i,j)]); end
    end
  end
  
  adjk=adjk*adj;
 
end

clear adjk

betw = zeros(1,n);
for i=1:n
    [dist,P]=dijkstra(adj,i,[]);

    for j=1:n
      
        if dist(j)<=1; continue; end   % either i=j or i,j are 1 edge apart
        betw(P{j}(2:dist(j))) = betw(P{j}(2:dist(j))) + 1/spaths(i,j);
        
    end
end

betw=betw/(2*nchoosek(n,2));   % further normalize by the number of all node pairs