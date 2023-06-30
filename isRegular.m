% Check whether a graph is regular, i.e. whether every node has the same degree.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean, 0 or 1
%
% Note: Defined for unweighted graphs only.
% Last updated: Sep 23, 2012

function S = isRegular(adj)

    S = false;

    degs = sum(adj > 0); % remove weights and sum columns

    if degs == degs(1) * ones(size(degs));
        S = true;
    end


%!test
%! adj = edgeL2adj(canonicalNets(20,'cycle'));
%! assert(isRegular(adj),true)

%! adj = edgeL2adj(canonicalNets(20,'tree',3));
%! assert(isRegular(adj),false)

%! assert(isRegular([0 1; 1 0]),true)
%! assert(isRegular([0 0; 1 0]),false)

%!demo
%! % undirected binary tree
%! adj = [0 1 1; 1 0 0; 1 0 0];
%! isRegular(adj)

%! % undirected 3−node cycle
%! adj = [0 1 1; 1 0 1; 1 1 0]; 
%! isRegular(adj)

%! % same as above, but edges are weighted
%! adj = [0 2 2; 2 0 2; 2 2 0];
%! isRegular(adj)

%! % a 4−node cycle
%! isRegular([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0])