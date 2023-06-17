% Extract the giant component of a graph;
% The giant component is the largest connected component.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: giant component matrix and node indices of the giant component
%
% Other routines used: findConnComp.m, subgraph.m
% Last updated: September 22, 2012

function [GC, gc_nodes] = giantComponent(adj)

    comp = findConnComp(adj);

    L = [];

     % computing component sizes
    for k = 1:length(comp);
        L = [L, length(comp{k})];
    end

    [maxL, ind_max] = max(L);

    gc_nodes = comp{ind_max};
    GC = subgraph(adj, gc_nodes);



%!test
%!assert(giantComponent([0 1 0; 1 0 0; 0 0 0]),[0 1; 1 0])
%! clear modules
%! modules{1}=[0];
%! randint = randi(10)+1;
%! Adj = []; adj = [];
%! % make up a matrix (Adj) of randint disconnected components (adj)
%! for x=1:randint
%!   randsecint = randi(10)+5;
%!   lastnode = modules{length(modules)}(length(modules{length(modules)}));
%!   modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
%!   % make sure adj is not empty, is connected and the number of nodes is "randsecint"
%!   while isempty(adj) || not(isConnected(adj)) || not(length(adj)==randsecint)
%!     adj=randomGraph(randsecint,0.5); 
%!   end
%!   Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
%! end
%! modules=modules(2:length(modules));
%! L = [];
%! for m=1:length(modules); 
%!   L = [L, length(modules{m})]; 
%! end;
%! [maxL,maxind] = max(L);
%! [GC, GCnodes] = giantComponent(Adj);
%! assert(GC, subgraph(Adj,modules{maxind}))
%! assert(GCnodes, modules{maxind})

%!demo
%! adj = [0 1 0; 1 0 0; 0 0 1];
%! [GC, I] = giantComponent(adj)
%! adj = [0 1 2 0; 1 0 0 0; 2 0 0 0; 0 0 0 0];
%! [GC, I] = giantComponent(adj)
