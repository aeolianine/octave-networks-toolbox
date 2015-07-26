% This function returns the betweenness measure of all vertices.
% Betweenness centrality measure: number of shortest paths running through a vertex.
% 
% Note 1: Valid for a general graph (multiple shortest paths possible).
% Note 2: Currently this function is quite slow, because of
%         findAllShortestPaths.m. This is especially true for dense
%         graphs.
%
% INPUTS: adjacency or distances matrix (nxn)
% OUTPUTS: betweeness vector for all vertices (1xn)
%
% Other routines used: numNodes.m, adj2adjL.m, simpleDijkstra.m, findAllShortestPaths.m
% GB: July 3, 2014

function betw = nodeBetweenness(adj)

n = numNodes(adj);
spaths=inf(n,n);
betw = zeros(1,n);  
adjL = adj2adjL(adj);

for i=1:n
    uB = simpleDijkstra(adj,i);
    for j=1:n
        if i==j; continue; end

        [allPaths, ~] = findAllShortestPaths(adjL,i,j, uB(j), allPaths={},path=[]);
        spaths(i,j) = length(allPaths);
        
        % for all paths in allPaths, parse out the path:
        for p=1:length(allPaths)
            path = strsplit(allPaths{p},'-');
            pathvec = [];
            for x=2:length(path)  % skip the first one
                pathvec = [pathvec str2num(path{x})];
            end
            betw(pathvec(2:length(pathvec)-1)) = betw(pathvec(2:length(pathvec)-1)) + 1/spaths(i,j);
        end
        
    end  % end of j=1:n
end      % end of i=1:n
    
% this last step is just additional normalization, and is arbitrary 
betw=betw/(2*nchoosek(n,2));