% Implementation of a community finding algorithm by Blondel et al
% Source: "Fast unfolding of communities in large networks", July 2008
%          https://sites.google.com/site/findcommunities/
% Note 1: This is just the first step of the Louvain community finding
%         algorithm. To extract fewer communities, need to repeat with 
%         the resulting modules themselves.
% Note 2: Works for undirected graphs only.
% Note 3: Permuting randomly the node order at every step helps the
%         algorithm performance. Unfortunately, node order in this
%         algorithm affects the results.
%
% INPUTs: adjancency matrix, nxn
% OUTPUTs: modules, and node community labels (inmodule)
%
% Other routines used: numEdges.m, numNodes.m, kneighbors.m
% Last updated, Oct 17 2012

function [modules, inmodule] = louvainCommunityFinding(adj)

    m = numEdges(adj);
    n = numNodes(adj);

    inmodule = {}; % inmodule{node} = module

    for mm = 1:length(adj)
        inmodule{mm} = mm;
    end

    modules = inmodule; % equal only for this step; modules{ind} = [nodes in module]

    % for all nodes, visit and assess joining to their neighbors
    % revisit until no improvement in dQ is available

    converged = false;

    while not(converged)

        converged = true;

        perm = randperm(n); % permute the order of the nodes randomly
        % this helps the algorithm performance
        for i = 1:n

            i = perm(i);

            neigh = kneighbors(adj, i, 1);
            dQ = zeros(1, length(neigh));

            c0 = inmodule{i};

            for nei = 1:length(neigh)
                % attempt to join i and neigh(nei)
                % that is: move i from modules{i} to modules{neigh(nei)}
                % if dQs balance is positive, do move; else: move on

                c1 = inmodule{neigh(nei)};
                dQ(nei) = dQi(i, c0, c1, adj, modules, m);

            end % loop across all neighbors

            [maxdQ, maxnei] = max(dQ);

            if maxdQ > 0

                % remove i from c0=inmodule{i}
                modules{inmodule{i}} = setdiff(modules{inmodule{i}}, i);
                % assign inmodule{i} to inmodule{neigh(maxnei)}
                inmodule{i} = inmodule{neigh(maxnei)};
                % add i to modules{inmodule{neigh(maxnei)}}
                modules{inmodule{neigh(maxnei)}} = [modules{inmodule{neigh(maxnei)}} i];

                converged = false;

            end

        end % loop across all nodes

    end % convergence loop

    % remove empty modules
    new_modules = {};

    for mm = 1:length(modules)

        if length(modules{mm}) > 0
            new_modules{length(new_modules) + 1} = modules{mm};
        end

    end

    modules = new_modules;

    printf('found %3i modules\n', length(modules));


function dqi = dQi(i, c0, c1, adj, modules, m)

    % INPUTs:
    %        c0, c1: indices of the two modules
    %        i: node which changes membership
    %        modules: modules{c0} = [list of nodes] etc
    %        m: total number of nodes in graph
    %
    % OUTPUTs: dqi: modularity change if i moves from c0 to c1

    if c0 == c1
        dqi = 0;
        return;
    end % same module - no change in modularity

    % k_ic0: degree of i within c0
    k_ic0 = sum(adj(i, modules{c0}));

    % k_ic1: degree of i witnin c1
    k_ic1 = sum(adj(i, modules{c1}));

    % k_i: degree of i (in entire graph)
    k_i = sum(adj(i, :));

    % m_c0: number of edges in c0
    m_c0 = numEdges(adj(modules{c0}, modules{c0}));

    % m_c1: number of edges in c1
    m_c1 = numEdges(adj(modules{c1}, modules{c1}));

    dqi = (k_ic1 - k_ic0) / (2 * m) - k_i * (m_c1 - m_c0) / (2 * m^2);


function p = randperm(n)
    [~, p] = sort(rand(n, 1));



%!test
%! extended_bowtie0 = [1 2; 2 3; 3 2; 3 4; 4 5; 5 6; 4 6; 6 7; 7 8; 7 9; 8 9];
%! extended_bowtie = [];
%! for row=1:size(extended_bowtie0,1)
%!   extended_bowtie = [extended_bowtie; extended_bowtie0(row,:) 1];
%! end
%! clear extended_bowtie0
%! extended_bowtie = symmetrizeEdgeL(extended_bowtie);
%! adj = edgeL2adj(extended_bowtie);

%! [modules,inmodule]=louvainCommunityFinding(adj);
%! assert(length(modules),3)
%! assert([inmodule{1},inmodule{2},inmodule{3}],[1,1,1]*inmodule{1})
%! assert([inmodule{4},inmodule{5},inmodule{6}],[1,1,1]*inmodule{4})
%! assert([inmodule{7},inmodule{8},inmodule{9}],[1,1,1]*inmodule{7})

%!test
%! T = load_test_graphs();
%! bowtie = T{4}{2};
%! [modules,inmodule]=louvainCommunityFinding(bowtie);
%! assert(length(modules),2)
%! assert([inmodule{1},inmodule{2},inmodule{3}],[1,1,1]*inmodule{1})
%! assert([inmodule{4},inmodule{5},inmodule{6}],[1,1,1]*inmodule{4})

%!test
%! % concatenate 4 complete graphs: 
%! adj = ones(10,10)-eye(10);
%! Adj = zeros(40,40);
%! 
%! Adj(1:10,1:10)=adj;
%! Adj(11:20,11:20)=adj;
%! Adj(21:30,21:30)=adj;
%! Adj(31:40,31:40)=adj;
%! 
%! Adj(10,11) = 1; Adj(11,10) = 1;
%! Adj(20,21) = 1; Adj(21,20) = 1;
%! Adj(30,31) = 1; Adj(31,30) = 1;
%! 
%! [modules,inmodule]=louvainCommunityFinding(Adj);
%! assert(length(modules),4)

%!test
%! % concatenate 4 dense graphs
%! adj = [0 1; 0 0];
%! while not(isConnected(adj)); adj = randomGraph(10,0.9); end
%! Adj = zeros(40,40);
%! 
%! Adj(1:10,1:10)=adj;
%! Adj(11:20,11:20)=adj;
%! Adj(21:30,21:30)=adj;
%! Adj(31:40,31:40)=adj;
%! 
%! Adj(10,11) = 1; Adj(11,10) = 1;
%! Adj(20,21) = 1; Adj(21,20) = 1;
%! Adj(30,31) = 1; Adj(31,30) = 1;
%! 
%! [modules,inmodule]=louvainCommunityFinding(Adj);
%! assert(length(modules),4)


%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! [modules,inmodule]=louvainCommunityFinding(bowtie);
%! assert(length(modules),2)
%! modules
%! inmodule