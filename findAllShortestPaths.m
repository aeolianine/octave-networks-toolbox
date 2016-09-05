% Return all shortest paths from a start to an end node, given a graph.
%
% INPUTS: adjacency list, adjL (1xn), start node "s" (index between 1 and n), 
%         end node "t" (index between 1 and n)
%         upperBound: max path length allowed in the search; best
%         to set this to the graph diameter; if many shortest paths
%         are sought at the same time, or the graph is very dense, 
%         upperBound can be set to the length of one pre-computed 
%         shortest path, using simpleDijkstra.m for example.
% OUTPUTS: list of all shortest paths between "s" and "t"
%
% Note 1: Works for a un/directed, unweighted graph.
% Note 2: This function uses recursion. It can be quite slow for
%         dense graphs.
% Other routines used: findAllShortestPaths.m
% GB: last updated, July 15, 2015

function [allPaths, upperBound] = findAllShortestPaths(adjL,s,t,upperBound, allPaths={}, path = [])


path = [path s];  % add "s" to path

if s == t
    if length(path)-1<=upperBound
        upperBound = length(path)-1;
        pathstr = '';
        for i=1:length(path); pathstr = strcat(pathstr,'-',num2str(path(i))); end
        allPaths{length(allPaths) + 1} = pathstr;
    end
else
    if length(path)-1<=upperBound
        for node=1:length(adjL{s})
            node = adjL{s}(node);
            if length(find(node==path))==0  % avoid cycles
                [allPaths,upperBound] = findAllShortestPaths(adjL, node, t, upperBound, allPaths, path);
            end
        end
    end
end

    
path2remove = {};
for i=1:length(allPaths)
    pathstr = allPaths{i};
    if length(strsplit(pathstr,'-'))-2>upperBound
        path2remove{length(path2remove)+1} = pathstr;
    end
end

allPaths = setdiff(allPaths,path2remove);