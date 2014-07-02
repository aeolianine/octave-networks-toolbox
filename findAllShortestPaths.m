% Return all shortest paths from a start to an end node, given a graph.
%
% Note: This function uses recursion. 
%
% INPUTS: adjacency list, adjL (1xn), start node "s" (index between 1 and n), 
%         end node "t" (index between 1 and n)
% OUTPUTS: list of all shortest paths between "s" and "t"
%
% Note: Works for a un/directed, unweighted graph.
% GB: last updated, June 30, 2014

function [allPaths, upperBound] = findAllShortestPaths(adjL,s,t, allPaths={}, path = [], upperBound = length(adjL)-1)


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
                [allPaths,upperBound] = findAllShortestPaths(adjL, node, t, allPaths, path, upperBound);
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