% Check whether a graph is a tree.
% A tree is a connected graph with n nodes and (n-1) edges.
% Source: "Intro to Graph Theory" by Bela Bollobas
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean variable, 0 or 1
%
% Other routines used: isConnected.m, numEdges.m, numNodes.m
% Last updated: Sep 24, 2012

function S = isTree(adj)

    S = false;

    if isConnected(adj) && numEdges(adj) == numNodes(adj) - 1;
        S = true;
    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(isTree(T{1}{2}), false)
%!assert(isTree(T{2}{2}), true)
%!assert(isTree(T{3}{2}), false)
%!assert(isTree(T{4}{2}), false)
%!assert(isTree(T{5}{2}), false)
%!assert(isTree(edgeL2adj(T{10}{2})), true)
%!assert(isTree(edgeL2adj(T{11}{2})), false)
%!assert(isTree(T{13}{2}), false)
%!assert(isTree(T{14}{2}), false)
%!assert(isTree(T{16}{2}), false)
%!assert(isTree(T{18}{2}), false)
%!assert(isTree(edgeL2adj(T{19}{2})), true)

%! adj = edgeL2adj(canonicalNets(randi(10)+10,'hexlattice'));
%! assert(isTree(adj),false)
%! adj = edgeL2adj(canonicalNets(randi(10)+10,'trilattice'));
%! assert(isTree(adj),false)
%! adj = edgeL2adj(canonicalNets(randi(10)+10,'hierarchy', b=3));
%! assert(isTree(adj),false)

%!demo
%! isTree([0 1 1; 1 0 0; 1 0 0])
%! isTree([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0])