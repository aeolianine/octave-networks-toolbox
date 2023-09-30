% Return the number of nodes, given an adjacency list, or adjacency matrix
% INPUTs: adjacency list: {i:j_1,j_2 ..} or adjacency matrix, ex: [0 1; 1 0]
% OUTPUTs: number of nodes, integer
%
% Last update: Sep 19, 2012

function n = numNodes(adjL)

    n = length(adjL);


%!test
%! randint = randi(101);
%! assert(numNodes(randomGraph(randint)),randint)

%!test
%!shared T
%! T = load_test_graphs();
%! for i=1:length(T)
%!     if strcmp(T{i}{3},'adjacency')
%!         assert( numNodes(T{i}{2}), T{i}{6} )
%!     end
%! end

%!demo
%! n = randi(100);
%! adj = randomGraph(n);
%! assert(numNodes(adj), n)
%! adjL = {[2, 3], [1, 3], [1, 2, 4], [3, 5, 6], [4, 6], [4, 5]};
%! assert(numNodes(adjL), 6)
