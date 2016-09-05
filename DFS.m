% Find all paths from a start to an end node, 
% bounded by a constant, using depth-first-search.
%
% Source: Idea from 6.002x, Spring 2014.
% Note: Uses recursion.
% 
% INPUTs: graph structure (a nxn adjacency matrix)
%         s - start node
%         e - end node
%         upperBound = 5, some constant bounding the 
%                    length of the path; typically size(adj,1)-1
%         allPaths = {}, path = []; are by default empty, 
%                                  serving the recursion
%
% OUTPUTs: a list {} of all shortest paths from "s" to "e", 
%                                shorter than "upperBound"
%
% Other routines used: kneighbors.m, DFS.m
% GB: last updated Oct 27 2014


function allPaths = DFS(graph, s, e, allPaths = {}, path = [], upperBound = 5)

        path = [path s];
        
        if s == e
            if length(path)-1<=upperBound   % this can be modified to include edge weights [(1,2,{}),..]
                allPaths{length(allPaths)+1} = path;
            end
        else
            if length(path)-1<=upperBound
              
                kneigh = kneighbors(graph, s, 1);
                
                for node=1:length(kneigh)
                    
                    node = kneigh(node);
                    
                    if length(find(path==node))==0   % avoid cycles
                        allPaths = DFS(graph, node, e, allPaths, path, upperBound);
                    end
                    
                end
            end
        end
