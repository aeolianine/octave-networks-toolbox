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
% Last updated: Oct 27 2014

function allPaths = DFS(graph, s, e, allPaths = {}, path = [], upperBound = 5)

    path = [path s];

    if s == e

        if length(path) - 1 <= upperBound % this can be modified to include edge weights [(1,2,{}),..]
            allPaths{length(allPaths) + 1} = path;
        end

    else

        if length(path) - 1 <= upperBound

            kneigh = kneighbors(graph, s, 1);

            for node = 1:length(kneigh)

                node = kneigh(node);

                if length(find(path == node)) == 0 % avoid cycles
                    allPaths = DFS(graph, node, e, allPaths, path, upperBound);
                end

            end

        end

    end



%!test
%!shared allPaths, T
%! T = load_test_graphs();

%! allPaths = DFS(T{1}{2}, 1, 2, allPaths = {}, path = [], upperBound = 2);
%!assert(allPaths, {[1 2]})
%! allPaths = DFS(T{1}{2}, 2, 1, allPaths = {}, path = [], upperBound = 2);
%!assert(allPaths, {})

%! allPaths = DFS(T{2}{2}, 1, 2, allPaths = {}, path = [], upperBound = 2);
%!assert(allPaths, {[1 2]})
%! allPaths = DFS(T{2}{2}, 2, 1, allPaths = {}, path = [], upperBound = 2);
%!assert(allPaths, {[2 1]})

%! allPaths = DFS(T{4}{2}, 1, 5, allPaths = {}, path = [], upperBound = 1);
%!assert(allPaths, {})
%! allPaths = DFS(T{4}{2}, 1, 5, allPaths = {}, path = [], upperBound = 3);
%!assert(allPaths, {[1 3 4 5]})

%! allPaths = DFS(T{5}{2}, 1, 5, allPaths = {}, path = [], upperBound = 4);
%!assert(allPaths, {})
%! allPaths = DFS(T{5}{2}, 1, 3, allPaths = {}, path = [], upperBound = 3);
%!assert(allPaths, {[1 2 3], [1 3]})

%! allPaths = DFS(T{5}{2}, 1, 3, allPaths = {}, path = [], upperBound = 1);
%!assert(allPaths, {[1 3]})

%! allPaths = DFS([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0], 1, 3, allPaths = {}, path = [], upperBound = 3);
%!assert(allPaths, {[1 2 3], [1 4 3]})
%! allPaths = DFS([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0], 4, 2, allPaths = {}, path = [], upperBound = 3);
%!assert(allPaths, {[4 1 2], [4 3 2]})


%!demo
%! adj = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
%! allPaths = DFS(adj , 1, 3, allPaths = {}, path = [] , upperBound = 3)
%! % with a smaller upperBound , the paths are not found
%! allPaths = DFS(adj , 1, 3, allPaths = {}, path = [] , upperBound = 1)