% Build a random modular graph, given number of modules, and link densities.
%
% INPUTs: number of nodes (n), number of modules (c), total link density (p),
%         and ratio of nodal degree to nodes within the same module
%         to the degree to nodes in other modules (r);
%         if specified, "labels" overwrites the random node
%         assignment to clusters, eg: labels = [1,1,2,3,3,4]
% OUTPUTs: adjacency matrix, modules to which the nodes are assigned
%
% Idea and code about pre-specified labels by Jonathan Hadida, June 12, 2014
% Last updated: July 6, 2014

function [adj, modules] = randomModularGraph(n, c, p, r, labels)

    % n - number of nodes
    % c - number of clusters/modules
    % p - overall probability of attachment
    % r - proportion of links within modules
    % labels - pre-specified cluster assignments

    % assign nodes to modules: 1 -> n/c, n/c+1 -> 2n/c, ... , (c-1)n/c -> c(n/c);

    if nargin < 5

        % equal, up to rounding assignment of nodes to clusters
        modules = {};

        for k = 1:c
            modules{k} = round((k - 1) * n / c + 1):round(k * n / c);
        end

    elseif nargin == 5

        assert(length(labels) == n)

        % pre-specified clusters with "labels"
        md = [0, find(diff(labels)), n]';
        md = [md(1:end - 1) + 1, md(2:end)];
        c = length(unique(labels)); % overwrite "c", the number of clusters, if needed

        modules = {};

        for k = 1:c
            modules{k} = md(k, 1):md(k, 2);
        end

    end

    adj = zeros(n); % initialize adjacency matrix

    % DERIVATION of probabilities
    % k_in/k_out = r, k_in + k_out = k = p(n-1)
    % => (1/r)k_in + k_in = p(n-1), => k_in = rp(n-1)/(r+1), k_out = p(n-1)/(r+1)
    % k_in = p_in*(n/c-1) => p_in = rpc(n-1)/((r+1)(n-c))
    % k_out = p_out*(n-n/c) => p_out = pc(n-1)/(n(r+1)(c-1))

    p_in = r * p * c * (n - 1) / ((r + 1) * (n - c));
    p_out = p * c * (n - 1) / (n * (r + 1) * (c - 1));

    for i = 1:n

        for j = i + 1:n

            module_i = ceil(c * i / n); % the module to which i belongs to
            module_j = ceil(c * j / n); % the module to which j belongs to

            if module_i == module_j

                % prob of attachment within module
                if rand < p_in
                    adj(i, j) = 1;
                    adj(j, i) = 1;
                end

            else

                % prob of attachment across modules
                if rand < p_out
                    adj(i, j) = 1;
                    adj(j, i) = 1;
                end

            end

        end

    end



%!test
%! for x=1:10
%!   N = randi(50)+10;
%!   c = randi(5)+1;
%!   [adj, modules] = randomModularGraph(N,c,0.2,4);
%!   assert(numNodes(adj),N)
%!   assert(length(modules),c)
%! end

%!test
%! % ..... testing with fixed labels ................
%! [adj, modules] = randomModularGraph(5,2,0.2,0.9,[1,1,2,3,4]);
%! mods = {}; mods{1} = [1,2]; mods{2} = [3]; mods{3} = [4]; mods{4} = 5;
%! assert(modules, mods)
%! [adj, modules] = randomModularGraph(6,0,0.2,0.8,[1,1,1,2,2,2]);
%! mods = {}; mods{1} = [1,2,3]; mods{2} = [4,5,6];
%! assert(modules, mods)
%! [adj, modules] = randomModularGraph(4,0,0.2,0.5,[1,2,3,4]);
%! mods = {}; mods{1} = 1; mods{2} = 2; mods{3} = 3; mods{4} = 4;
%! assert(modules, mods)
%! [adj, modules] = randomModularGraph(4,0,0.2,0.5,[1,1,1,1]);
%! mods = {}; mods{1} = [1,2,3,4];
%! assert(modules, mods)
%! [adj, modules] = randomModularGraph(75,4,0.2,0.5,[ones(1,20),2*ones(1,20),3*ones(1,20),4*ones(1,15)]);
%! mods = {}; mods{1} = [1:20]; mods{2} = [21:40]; mods{3} = [41:60]; mods{4} = [61:75];
%! assert(modules, mods)


%!test
%! % ... testing the in/out-degree ratio distribution ......
%! ratio = [];
%! for x=1:100
%!   
%!   mods = 4;
%!   N = 100;
%!   dens = log(N)/N;  % threshold of connectivity
%!   
%!   adj = [];
%!   while not(isConnected(adj)); 
%!     [adj,modules] = randomModularGraph(N,mods,dens,10); 
%!   end
%!   
%!   for M=1:mods
%! 
%!     kin = [];
%!     kout = [];
%!     
%!     for node=1:length(modules{M})
%!       i = modules{M}(node);
%! 
%!       ss = 0;
%!       for c=1:length(modules)
% !        if c == M; continue; end
%!         ss += sum(adj(i,modules{c}));
%!       end
%!   
%!       if ss==0; ss=0.01; end
%!       
%!       kin = [kin sum(adj(i,modules{M}))];
%!       kout = [kout ss];
%! 
%!       
%!     end
%!     
%!     kin = sum(kin)/length(kin);
%!     kout = sum(kout)/length(kout);
%!     
%!     ratio = [ratio kin/kout];
%!   end
%!  
%!   assert(length(modules),mods)
%!   assert(length(adj),N)
%! end
%! 
%! hist(ratio,20)
%! title('supposed to be centered around 10')
%! hold off;


%!demo
%! [adj, modules] = randomModularGraph(100, 4, 0.1, 0.9);
%! assert(numNodes(adj), 100)
%! assert(length(modules), 4)
%! [adj, modules] = randomModularGraph(4, 4, 0.1, 0.5, [1, 1, 2, 2]);
%! assert(length(modules), 2)