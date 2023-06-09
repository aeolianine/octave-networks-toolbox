% Algorithm for finding connected components in a graph
% Note: Valid for undirected graphs only.
%
% INPUTS: adj - adjacency matrix, nxn
% OUTPUTS: a list of the components comp{i}=[j1,j2,...jk]
%
% Other routines used: findConnCompI.m, degrees.m
% Last updated: September 22, 2012

function comp_mat = findConnComp(adj)

    [deg, ~, ~] = degrees(adj); % degrees
    comp_mat = {}; % initialize components matrix

    for i = 1:length(deg)

        if deg(i) > 0
            done = 0;

            for x = 1:length(comp_mat)

                if length(find(comp_mat{x} == i)) > 0 % i in comp_mat(x).mat
                    done = 1;
                    break
                end

            end

            if not(done)
                comp = findConnCompI(adj, i);
                comp_mat{length(comp_mat) + 1} = comp;
            end

        elseif deg(i) == 0
            comp_mat{length(comp_mat) + 1} = [i];
        end

    end


%!test
%!shared T
%! T = load_test_graphs();
%!assert(findConnComp(T{5}{2}),{[1,2,3],[4,5,6]})

%! clear modules
%! modules{1}=[0];
%! randint = randi(21);
%! Adj = []; adj = [];
%! % make up a matrix (Adj) of randint disconnected components (adj)
%! for x=1:randint
%!   randsecint = randi(5)+5;
%!   
%!   % remember the disconnected components in "modules"
%!   lastnode = modules{length(modules)}(length(modules{length(modules)}));
%!   modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
%!   
%!   % make sure adj is not empty, is connected and the number of nodes is "randsecint"
%!   while isempty(adj) || not(isConnected(adj)) || not(length(adj)==randsecint); adj=randomGraph(randsecint,0.5); end
%! 
%!   Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
%! end
%! 
%! modules=modules(2:length(modules));
%! assert(findConnComp(Adj),modules)

%!demo
%! % two disconnected three−node cycles
%! adj=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! comp = findConnComp(adj)
%! adj = [0 1 1; 1 0 0; 1 0 0];
%! findConnComp(adj)