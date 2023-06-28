% Convert an incidence matrix to an edge list.
%
% Inputs: inc - incidence matrix nxm (number of nodes x number of edges)
% Outputs: edge list - mx3, m x (node 1, node 2, edge weight)
%
% Example: [-1; 1] <=> [1,2,1], one directed (1->2) edge
% Last updated: Sep 25 2012

function el = inc2edgeL(inc)

    m = size(inc, 2); % number of edges
    el = zeros(m, 3); % initialize edge list [n1, n2, weight]

    for e = 1:m
        ind_m1 = find(inc(:, e) == -1);
        ind_p1 = find(inc(:, e) == 1);

        if numel(ind_m1) == 0 && numel(ind_p1) == 1 % undirected, self-loop
            el(e, :) = [ind_p1 ind_p1 1];

        elseif numel(ind_m1) == 0 && numel(ind_p1) == 2 % undirected
            el(e, :) = [ind_p1(1) ind_p1(2) 1];
            el = [el; ind_p1(2) ind_p1(1) 1];

        elseif numel(ind_m1) == 1 && numel(ind_p1) == 1 % directed
            el(e, :) = [ind_m1 ind_p1 1];

        end

    end


%!test
%!assert(inc2edgeL([1 0 0; 0 1 0; 0 0 1]),[1 1 1; 2 2 1; 3 3 1])  % three self-loops
%!assert(inc2edgeL([-1 -1; 1 0; 0 1]),[1 2 1; 1 3 1])   % tree 3 nodes
%!assert(inc2edgeL([1 1; 1 1]),[1 2 1; 1 2 1; 2 1 1; 2 1 1])     % one double edge

%!shared T
%! T = load_test_graphs();
%!assert(inc2edgeL([-1;1]),T{1}{5})                     % one directed edge
%!assert(inc2edgeL([ 1;1]),T{2}{5})                     % one undirected edge
%!assert(sortrows(inc2edgeL(T{15}{2})(1:length(T{15}{5}),1:2)),sortrows(T{15}{5}))

%!demo
%! inc = [-1; 1]; % single directed edge
%! inc2edgeL(inc) 
%! inc = [1 0; 1 1; 0 1]; % two undirected edges
%! inc2edgeL(inc) 
%! inc = [-1 0; 1 -1; 0 1];  % two directed edges
%! inc2edgeL(inc) 
