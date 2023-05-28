% Convert adjacency matrix to an incidence matrix
% Note: Valid for directed/undirected, simple/not simple graphs
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: incidence matrix: n x m (number of edges)
%
% Other routines used: isDirected.m
% Last updated: May 20 2023

function inc = adj2inc(adj)

    n = length(adj); % number of nodes

    inc = []; % initialize incidence matrix

    if isDirected(adj)

        for i = 1:n

            for j = 1:n

                % handle self-loops
                if i == j; for x = 1:adj(i, j); inc = [inc; zeros(1, length(1:i - 1)), 1, zeros(1, length(i + 1:n))]; end; continue; end

                for x = 1:adj(i, j)% add multiple edges if any

                    if i < j
                        inc = [inc; zeros(1, length(1:i - 1)), -1, zeros(1, length(i + 1:j - 1)), 1, zeros(1, length(j + 1:n))];
                    else
                        inc = [inc; zeros(1, length(1:j - 1)), 1, zeros(1, length(j + 1:i - 1)), -1, zeros(1, length(i + 1:n))];
                    end

                end

            end

        end

    else % undirected

        for i = 1:n

            for j = i:n

                % handle self-loops
                if i == j; for x = 1:adj(i, j); inc = [inc; zeros(1, length(1:i - 1)), 1, zeros(1, length(i + 1:n))]; end; continue; end
                % add multiple edges if any
                for x = 1:adj(i, j); inc = [inc; zeros(1, length(1:i - 1)), 1, zeros(1, length(i + 1:j - 1)), 1, zeros(1, length(j + 1:n))]; end

            end

        end

    end

    inc = inc';


%!test
%!shared T, randint
%! T = load_test_graphs();
%! randint = randi(10)+1;
%!assert(adj2inc(eye(randint)),eye(randint))
%!assert(adj2inc(T{13}{2}), T{15}{2})   % directed 3-cycle
%!assert(adj2inc([0 1 0; 0 1 0; 1 0 0 ]),[-1 0 1; 1 1 0; 0 0 -1])
%!assert(adj2inc([0 2; 0 0]),[-1 -1; 1 1])  % directed double edge
%!assert(adj2inc(T{1}{2}), [-1 1]')          % one directed edge

%!demo
%! adj2inc([0 1 0; 1 0 0; 0 0 0])  % one undirected edge
%! adj2inc([0 1 1; 1 0 0; 1 0 0])  % two undirected edges
%! adj2inc([0 1 1; 0 0 0; 0 0 0])  % two directed edges
