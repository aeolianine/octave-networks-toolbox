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
% Last updated: July 15, 2015

function [allPaths, upperBound] = findAllShortestPaths(adjL, s, t, upperBound, allPaths = {}, path = [])

    path = [path s]; % add "s" to path

    if s == t

        if length(path) - 1 <= upperBound
            upperBound = length(path) - 1;
            pathstr = '';
            for i = 1:length(path); pathstr = strcat(pathstr, '-', num2str(path(i))); end
            allPaths{length(allPaths) + 1} = pathstr;
        end

    else

        if length(path) - 1 <= upperBound

            for node = 1:length(adjL{s})
                node = adjL{s}(node);

                if length(find(node == path)) == 0 % avoid cycles
                    [allPaths, upperBound] = findAllShortestPaths(adjL, node, t, upperBound, allPaths, path);
                end

            end

        end

    end

    path2remove = {};

    for i = 1:length(allPaths)
        pathstr = allPaths{i};

        if length(strsplit(pathstr, '-')) - 2 > upperBound
            path2remove{length(path2remove) + 1} = pathstr;
        end

    end

    allPaths = setdiff(allPaths, path2remove);


%!test
%!shared adjL, shortestPathLength, allPaths
%! adjL = {}; adjL{1} = [2]; adjL{2} = [];  % 1-2 edge
%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,1,2, 2, allPaths={});
%!assert(shortestPathLength,1)
%!assert(allPaths{1},'-1-2')
%!assert(length(allPaths),1)

%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,2,1, 2, allPaths={});
%!assert(shortestPathLength,2)  % not updated, equal to length(adjL)-1
%!assert(allPaths, {})
%!assert(length(allPaths),0)

%! % two alternative paths (1-2-4 and 1-3-4)
%! adjL = {}; adjL{1} = [2,3]; adjL{2} = [4]; adjL{3} = [4]; adjL{4}=[];
%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,1,4, 5, allPaths={});
%!assert(shortestPathLength,2)
%!assert(length(allPaths),2)
%!assert(allPaths{1},'-1-2-4')
%!assert(allPaths{2},'-1-3-4')

%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,4,2, 4, allPaths={});
%!assert(shortestPathLength,4) % not updated, equal to length(adjL)-1
%!assert(length(allPaths),0)

%! % a one-directional cycle
%! adjL={}; adjL{1}=[2]; adjL{2}=[3]; adjL{3}=[4]; adjL{4}=[1];
%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,1,4, 4, allPaths={});
%!assert(shortestPathLength,3)
%!assert(length(allPaths),1)
%!assert(allPaths{1},'-1-2-3-4')

%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,4,3, 4, allPaths={});
%!assert(shortestPathLength,3)
%!assert(length(allPaths),1)
%!assert(allPaths{1},'-4-1-2-3')

%! % undirected cycle
%! adjL={}; adjL{1}=[2,4]; adjL{2}=[3,1]; adjL{3}=[2,4]; adjL{4}=[1,3];
%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,2,4, 4, allPaths={});
%!assert(shortestPathLength,2)
%!assert(length(allPaths),2)
%!assert(allPaths{1},'-2-1-4')
%!assert(allPaths{2},'-2-3-4')

%! [allPaths, shortestPathLength] = findAllShortestPaths(adjL,3,1, 3, allPaths={});
%!assert(shortestPathLength,2)
%!assert(length(allPaths),2)
%!assert(allPaths{1},'-3-2-1')
%!assert(allPaths{2},'-3-4-1')

%!demo
%! L = {[2, 3], [1, 3], [1, 2, 4], [3, 5, 6], [4, 6], [4, 5]}; 
%! findAllShortestPaths(L, 1, 6, 3)
%! findAllShortestPaths(L, 5, 2, 3)
%! L = {[2,4], [1,3], [2,4], [3,1]};
%! findAllShortestPaths(L, 1, 3, 3)
%! findAllShortestPaths(L, 4, 2, 3)