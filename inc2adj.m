% Convert an incidence matrix representation to an adjacency matrix representation for an arbitrary graph.
%
% INPUTs: incidence matrix, nxm (num nodes x num edges)
% OUTPUTs: adjacency matrix, nxn
%
% Last updated: Jul 27, 2014

function adj = inc2adj(inc)

    m = size(inc, 2); % number of edges
    adj = zeros(size(inc, 1)); % initialize adjacency matrix

    if isempty(find(inc == -1)) % undirected graph

        for e = 1:m

            ind = find(inc(:, e) == 1);

            if length(ind) == 2
                adj(ind(1), ind(2)) = 1;
                adj(ind(2), ind(1)) = 1;
            elseif length(ind) == 1% selfloop
                adj(ind, ind) = 1;
            else
                fprintf('inc2adj(): invalid incidence matrix.\n')
                return
            end

        end

    else % directed graph (there are "-1"
        % entries in the incidence matrix)

        for e = 1:m

            ind1 = find(inc(:, e) == 1);
            indm1 = find(inc(:, e) == -1);

            if isempty(indm1) && length(ind1) == 1 % selfloop
                adj(ind1, ind1) = adj(ind1, ind1) + 1;
            elseif length(indm1) == 1 && length(ind1) == 1
                adj(indm1, ind1) = adj(indm1, ind1) + 1;
            else
                fprintf('inc2adj(): invalid incidence matrix.\n')
                return
            end

        end

    end


%!test
%!shared T
%! T = load_test_graphs();
%! randint = randi(10)+1;
%! assert(inc2adj(eye(randint))==eye(randint))

%!assert(inc2adj(ones(3) - eye(3)),ones(3) - eye(3))

%! % 1->2, 2->2, 3->1
%!assert([0 1 0; 0 1 0; 1 0 0 ],inc2adj([-1 0 1; 1 1 0; 0 0 -1]))
%!assert([0 2; 0 0],inc2adj([-1 -1; 1 1]))  % directed double edge
%!assert(T{1}{2}, inc2adj([-1 1]'))          % one directed edge

%! % two edges (1->2, 3->1)
%!assert(inc2adj([-1 1; 1 0; 0 -1])==[0 1 0; 0 0 0; 1 0 0])

%!demo
%! % an example in which the incidence happens to equal the adjacency 
%! inc = [1 0 1; 1 1 0; 0 1 1];
%! inc2adj(inc)

%! % two directed edges
%! inc = [-1 -1; 1 0; 0 1];
%! inc2adj(inc)

%! % one directed edge
%! inc = [-1; 1]
%! inc2adj(inc)
