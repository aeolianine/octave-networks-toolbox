% Simple implementation of breadth-first-search of a graph.
% Returns a breadth-first-search tree.
%
% INPUTs: adjacency list (nxn), "adjL"
%         start node index, "s"
%         end node index, "t"
% OUTPUTs: BFS tree, in adjacency list {} format (directed)
%          (This is the tree that the algorithm walks in
%          search of the target. If the target is not found,
%          this tree is effectively a spanning tree of the
%          entire graph.
%
% Last updated: Nov 8 2014

function T = BFS(adjL, s, t)

    discovered = [s]; % mark visited
    q = [s]; % add to queue
    T = cell(1, length(adjL)); % initialize search path/tree

    if s == t
        printf('BFS(): The start and the end node are the same.\n')
        return
    end

    while not(isempty(q))
        j = q(1); q = q(2:length(q)); % pop the front
        neigh = adjL{j};

        for nn = 1:length(neigh)

            if isempty(find(discovered == neigh(nn)))
                T{j} = [T{j}, neigh(nn)];

                if neigh(nn) == t
                    return
                end

                discovered = [discovered, neigh(nn)];
                q = [q, neigh(nn)];
            end

        end

    end

    printf('BFS(): Target node not found.\n')



%!test
%!shared adjL, T, tt
%! T = load_test_graphs();
%! adjL = {1:2, 2:[]};
%! tt = BFS(adjL, 1, 1);
%!assert(tt, {[], []})

%! tt = BFS(adjL, 1, 2);
%!assert(tt{1}, 2)
%!assert(length(tt),2)
%!assert(class(tt),'cell')

%! tt = BFS(adjL, 2, 2);
%!assert(tt, {[],[]})

%! tt = BFS(adjL, 2, 3);
%!assert(tt, {[],[]})
%!assert(length(tt),2)
%!assert(class(tt),'cell')

%! tt = BFS(adjL, 1, 3);
%!assert(tt{1}, 2)
%!assert(tt{2}, [])
%!assert(length(tt),2)
%!assert(class(tt),'cell')

%! tt = BFS(T{9}{2}, 1, 4);
%!assert(tt{1},[2 3])
%!assert(tt{2},[])
%!assert(tt{3},[4])
%!assert(tt{4},[])
%!assert(tt{5},[])
%!assert(tt{6},[])

%! tt = BFS(T{9}{2}, 2, 6);
%!assert(tt{2},[1 3])
%!assert(tt{1},[])
%!assert(tt{3},[4])
%!assert(tt{4},[5 6])
%!assert(tt{5},[])
%!assert(tt{6},[])

%! tt = BFS(T{9}{2}, 5, 2);
%!assert(tt{5},[4 6])
%!assert(tt{6},[])
%!assert(tt{4},[3])
%!assert(tt{3},[1 2])
%!assert(tt{2},[])
%!assert(tt{1},[])

%! tt = BFS(T{9}{2}, 5, 10);
%!assert(tt{5},[4 6])
%!assert(tt{6},[])
%!assert(tt{4},[3])
%!assert(tt{3},[1 2])
%!assert(tt{2},[])
%!assert(tt{1},[])


%!demo
%! % adjacency list representation of the bowtie graph (I>−<I)
%! L = {[2, 3], [1, 3], [1, 2, 4], [3, 5, 6], [4, 6], [4, 5]}; 
%! BFS(L, 1, 6)
%! BFS(L, 5, 3)