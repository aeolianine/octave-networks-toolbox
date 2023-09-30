% Computing the modularity for a given module/commnunity partition.
% Defined as: Q=sum_over_modules_i (eii-ai^2) (eq 5) in Newman and Girvan.
% eij = fraction of edges that connect community i to community j, ai=sum_j (eij)
%
% Source: Newman, Girvan, "Finding and evaluating community structure in networks"
%         Newman, "Fast algorithm for detecting community structure in networks"
%
% INPUTs: adjacency matrix, nxn
%         set of modules as cell array of vectors, ex: {[1,2,3],[4,5,6]}
% OUTPUTs: modularity metric, in [-1,1]
%
% Note: This computation makes sense for undirected graphs only.
% Other functions used: numEdges.m
% Last updated: October 16, 2012

function Q = modularityMetric(modules, adj)

    nedges = numEdges(adj); % total number of edges

    Q = 0;

    for m = 1:length(modules)

        e_mm = numEdges(adj(modules{m}, modules{m})) / nedges;
        a_m = sum(sum(adj(:, modules{m}))) / nedges - e_mm;
        Q = Q + (e_mm - a_m^2);

    end


%  % alternative: Q = sum_ij { 1/2m [Aij-kikj/2m]delta(ci,cj) } = 
%  % = sum_ij Aij/2m delta(ci,cj) - sum_ij kikj/4m^2 delta(ci,cj) = 
%  % = sum_modules e_cc - sum_modules (kikj/4m^2) =
%  % = sum_modules (e_cc - ((sum_i ki)/2m)^2)
%  
%  n = numNodes(adj);
%  m = numEdges(adj);
%  
%  % define the inverse of modules: node "i" <- module "c" if "i" in module "c"
%  mod={};
%  for mm=1:length(modules)
%    for ii=1:length(modules{mm})
%      mod{modules{mm}(ii)}=modules{mm};
%    end
%  end
%  
%  Q = 0;
%  
%  for i=1:n
%    for j=1:n
%      
%     if not(isequal(mod(i),mod(j))); continue; end
%  
%     Q = Q + (adj(i,j) - sum(adj(i,:))*sum(adj(j,:))/(2*m))/(2*m);
%     
%    end
%  end


%!test
%! for i = 1:10
%! 
%!     adj = [0 1; 0 0];
%!     num_modules = randi([2, 5]);
%! 
%!     while not(isConnected(adj));
%!         adj = randomModularGraph(30, num_modules, 0.1, 3);
%!     end
%! 
%!     % compare to newmanCommFast
%!     [mH, Q1] = newmanCommFast(adj);
%!     close all;
%!     Q2 = [];
%! 
%!     for m = 1:length(mH);
%!         Q2 = [Q2 modularityMetric(mH{m}, adj)];
%!     end
%! 
%!     assert(Q1, Q2)
%! 
%!     % compare to the newman-girvan routine
%!     [modules0, ~, Q0] = newmanGirvan(adj, num_modules);
%!     assert(Q0, modularityMetric(modules0, adj))
%! 
%! end


%!demo
%! bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! Q = modularityMetric({[1, 2, 3], [4, 5, 6]} , bowtie)