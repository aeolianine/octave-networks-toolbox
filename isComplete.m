% Check whether a graph is complete, i.e. whether every node is linked to every other node.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean variable, true/false
%
% Note: Only defined for unweighted graphs.
% Last updated: Sep 23, 2012

function S = isComplete(adj)

    S = false; % default

    adj = adj > 0; % remove weights
    n = length(adj);

    % all degrees "n-1" or "n" or w/ n selfloops
    if sum(adj) == ones(1, n) * (n - 1) || sum(adj) == ones(1, n) * n;
        S = true;
    end


%!test
%!assert(isComplete([0 1; 1 0]),true)
%!shared T
%! T = load_test_graphs();
%!assert(isComplete(T{2}{2}),true)
%!assert(isComplete(T{3}{2}),true)
%!assert(isComplete(T{4}{2}),false)

%! randint = randi(10)+10;
%! adj = ones(randint)-eye(randint);
%! assert(isComplete(adj),true)

%!demo
%! isComplete([0 1 1; 1 0 0; 1 0 0])
%! isComplete([0 1 1; 1 0 1; 1 1 0])
