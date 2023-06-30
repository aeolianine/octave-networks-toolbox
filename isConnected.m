% Determine if a graph is connected
% Idea by Ed Scheinerman, circa 2006,
%         source: http://www.ams.jhu.edu/~ers/matgraph/
%         routine: matgraph/@graph/isconnected.m
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean variable, 0 or 1
%
% Note: This function works only for undirected graphs.
% Last updated: Sep 23 2012

function S = isConnected(adj)

    if not(isempty(find(sum(adj) == 0)));
        S = false;
        return;
    end

    n = length(adj);
    x = [1; zeros(n - 1, 1)]; % [1,0,...0] nx1 vector

    while 1
        y = x;
        x = adj * x + x;
        x = x > 0;

        if x == y
            break;
        end

    end

    S = true;

    if sum(x) < n
        S = false;
    end


% Alternative 1 ..........................................................
% If the algebraic connectivity is > 0 then the graph is connected
% a=algebraic_connectivity(adj);
% S = false; if a>0; S = true; end

% Alternative 2 ..........................................................
% Uses the fact that multiplying the adj matrix to itself k times give the
% number of ways to get from i to j in k steps. If the end of the
% multiplication in the sum of all matrices there are 0 entries then the
% graph is disconnected. Computationally intensive, but can be sped up by
% the fact that in practice the diameter is very short compared to n, so it
% will terminate in order of log(n)? steps.
% function S=isconnected(el):
%
%     S=false;
%
%     adj=edgeL2adj(el);
%     n=numnodes(adj); % number of nodes
%     adjn=zeros(n);
%
%     adji=adj;
%     for i=1:n
%         adjn=adjn+adji;
%         adji=adji*adj;
%
%         if length(find(adjn==0))==0
%             S=true;
%             return
%         end
%     end

% Alternative 3 ............................................................
% Find all connected components, if their number is 1, the graph is
% connected. Use find_conn_comp(adj).

%!test
%!shared T
%! T = load_test_graphs();
%!assert(isConnected(T{1}{2}),false)
%!assert(isConnected(T{2}{2}),true)
%!assert(isConnected(T{4}{2}),true)
%!assert(isConnected(T{5}{2}),false)
%!assert(isConnected(transpose(T{5}{2})),false)
%!assert(isConnected(T{13}{2}),true)
%!assert(isConnected(T{14}{2}),true)

%! for x=1:100
%!   % test the connected case
%!   adj = [0 1; 1 0];                % initialize
%!   N = randi(100)+5;
%!   while isConnected(adj) || sum(adj)==0
%!      adj = randomGraph(N,log(N)/N); 
%!   end
%!   assert(isConnected(adj),false)
%!   assert(abs(algebraicConnectivity(adj)) <10^(-10),true)
%!   assert(transpose(isConnected(adj)),false)
%!   assert(length(findConnComp(adj))>1, true)
  
%!   % test the not-connected case
%!   adj = [0 1; 0 0];                % initialize
%!   while not(isConnected(adj))
%!      adj = randomGraph(N,log(N)/N); 
%!   end
%!   assert(isConnected(adj),true)
%!   assert(algebraicConnectivity(adj)>0,true)
%!   assert(transpose(isConnected(adj)),true)
%!   assert(length(findConnComp(adj))==1, true)
%!   
%! end

%!demo
%! % undirected binary tree with 3 nodes 
%! isConnected([0 1 1; 1 0 0; 1 0 0])
%! % two disconnected 3−node cycles
%! adj=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0]; 
%! isConnected(adj)