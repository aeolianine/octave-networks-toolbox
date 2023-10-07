% Find the stronly connected components in a directed graph
% Source: Tarjan, R. E. (1972), "Depth-first search and linear graph algorithms",
%                                       SIAM Journal on Computing 1 (2): 146-160
% Wikipedia description: http://en.wikipedia.org/wiki/Tarjan's_strongly_connected_components_algorithm
%
% Input: graph, set of nodes and edges, in adjacency list format,
%        example: L{1}=[2], L{2]=[1] is a single (1,2) edge
% Outputs: set of strongly connected components, in cell array format
%
% Other routines used: strongConnComp.m
% Last updated: Sep 22, 2012

function [GSCC, v] = tarjan(L)

    GSCC = {};
    ind = 1; % node number counter
    S = []; % An empty stack of nodes

    for ll = 1:length(L)
        v(ll).index = [];
        v(ll).lowlink = [];
    end % initialize indices

    for vi = 1:length(L)

        if isempty(v(vi).index)
            [GSCC, S, ind, v] = strongConnComp(vi, S, ind, v, L, GSCC); % visit new nodes only
        end

    end


%!test
%! L = {}; L{1} = 2; L{2} = 1;
%! GSCC = tarjan(L);
%! assert(length(GSCC),1)
%! assert(GSCC{1},[1,2])

%!test
%! L = {}; L{1} = 2; L{2} = [];
%! GSCC = tarjan(L);
%! assert(length(GSCC),2)
%! assert(GSCC{1},[2])
%! assert(GSCC{2},[1])

%!test
%! L={}; L{1}=[2,3]; L{2}=[1]; L{3}=[1]; L{4}=[1]; % cherry tree (binary) + extra node
%! GSCC = tarjan(L);
%! assert(length(GSCC),2)
%! assert(GSCC{1},[1,2,3])
%! assert(GSCC{2},4)

%!test
%! L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2]; L{4}=[1]; % triangle with extra node
%! GSCC = tarjan(L);
%! assert(length(GSCC),2)
%! assert(GSCC{1},[1,2,3])
%! assert(GSCC{2},4)

%!test
%! L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2,4]; L{4}=[5,6]; L{5}=[4,6]; L{6}=[4,5];
%! GSCC = tarjan(L);
%! assert(length(GSCC),2)
%! assert(length(GSCC{1}),3)
%! assert(length(GSCC{2}),3)

%!test
%! L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2]; L{4}=[5,6]; L{5}=[4,6]; L{6}=[4,5];
%! GSCC = tarjan(L);
%! assert(length(GSCC),2)
%! assert(length(GSCC{1}),3)
%! assert(length(GSCC{2}),3)

%!test
%! L={}; L{2}=[1]; L{3}=[1]; L{4}=[1];
%! GSCC = tarjan(L);
%! assert(length(GSCC),4)

%!test
%! for iter=1:100  % random graph testing ....
%!   % testing undirected graphs
%!   adj = [0 1; 0 0];  % initialize so that the while loop does not break
%!   while not(isConnected(adj)); adj = randomGraph(randi(50)+1,rand); end
%! 
%!   % test that the GSCC contains one component, namely the entire graph
%!   L=adj2adjL(adj);
%!   GSCC = tarjan(L);
%!   assert(length(GSCC),1)
%!   assert(GSCC{1},[1:length(adj)])
%!   
%!   % testing directed graph
%!   adj=randomDirectedGraph(randi(50)+1,rand*0.1);
%!   L=adj2adjL(adj);
%!   GSCC = tarjan(L);
%!   
%!   
%!   if isConnected(adj) && isConnected(transpose(adj)) && length(adj)>0
%!     
%!     % there should be one component containing all nodes
%!     assert(length(GSCC),1)
%!     assert(GSCC{1},[1:length(adj)])
%!     
%!     
%!   else  % disconnected directed graph
%!     
%!     ll=[];
%!     for gg=1:length(GSCC); ll=[ll length(GSCC{gg})]; end;
%!     [ml,maxll]=max(ll);
%!     
%!     % the largest strongly connected component either is trivial (one node), or is connected
%!     assert(isConnected(adj(GSCC{maxll},GSCC{maxll})) | length(GSCC{maxll})==1)
%!     
%!     % adding any other node to the largest strongly connected component should not make the subgraph connected
%!     for ii=1:length(adj)
%!       if isempty(find(GSCC{maxll}==ii))
%!         
%!         tryGC = [GSCC{maxll}, ii];
%!         assert( not(isConnected(adj(tryGC,tryGC))) )
%!         
%!       end
%!       
%!     end
%!     
%!   end
%!   
%! end
