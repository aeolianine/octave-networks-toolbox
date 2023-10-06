% Test code for "Octave tools for Network Analysis"

clear all
close all


% Set of test graphs, in various formats
T = load_test_graphs();

% add another graph test here ....
% ................................................

% ................................................
% ... graph representation functions .............
% ................................................


% Testing symmetrize.m ...........................
printf('testing symmetrize.m\n')
tic
for i=1:20
  adj = randomDirectedGraph(randi(10)+3,rand);
  assert(isSymmetric(symmetrize(adj)),true)
end
assert(symmetrize(T{1}{2}),T{2}{2})
assert(symmetrize(edgeL2adj(T{11}{2})),edgeL2adj(T{10}{2}))
assert(symmetrize(T{16}{2}),T{13}{2})
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing symmetrizeEdgeL.m ......................
printf('testing symmetrizeEdgeL.m\n')
tic
for x=1:20
  adj = randomDirectedGraph(randi(20)+2,rand); % create a random adjacency
  el = adj2edgeL(adj);
  if isempty(el); continue; end
  elsym = symmetrizeEdgeL(el);
  adjsym = edgeL2adj(elsym);
  assert(isSymmetric(adjsym),true)
end

assert(sortrows(symmetrizeEdgeL(T{1}{5}))(1:2,1:2), sortrows(T{2}{5})(1:2,1:2))
assert(sortrows(symmetrizeEdgeL(T{6}{5}))(1:14,1:2), sortrows(T{4}{5})(1:14,1:2) )
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% ................................................
% ... basic network theory functions .............
% ................................................


% Testing selfLoops.m ...........................
printf('testing selfLoops.m\n')
tic
assert( selfLoops( edgeL2adj( T{8}{2} ) ), 1 )
assert( selfLoops( T{14}{2} ), 2 )
assert(selfLoops(T{4}{2}),0)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing tarjan.m ..............................
printf('testing tarjan.m\n')
tic
L = {}; L{1} = 2; L{2} = 1;
GSCC = tarjan(L);
assert(length(GSCC),1)
assert(GSCC{1},[1,2])

L = {}; L{1} = 2; L{2} = [];
GSCC = tarjan(L);
assert(length(GSCC),2)
assert(GSCC{1},[2])
assert(GSCC{2},[1])


L={}; L{1}=[2,3]; L{2}=[1]; L{3}=[1]; L{4}=[1]; % cherry tree (binary) + extra node
GSCC = tarjan(L);
assert(length(GSCC),2)
assert(GSCC{1},[1,2,3])
assert(GSCC{2},4)


L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2]; L{4}=[1]; % triangle with extra node
GSCC = tarjan(L);
assert(length(GSCC),2)
assert(GSCC{1},[1,2,3])
assert(GSCC{2},4)


L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2,4]; L{4}=[5,6]; L{5}=[4,6]; L{6}=[4,5];
GSCC = tarjan(L);
assert(length(GSCC),2)
assert(length(GSCC{1}),3)
assert(length(GSCC{2}),3)

L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2]; L{4}=[5,6]; L{5}=[4,6]; L{6}=[4,5];
GSCC = tarjan(L);
assert(length(GSCC),2)
assert(length(GSCC{1}),3)
assert(length(GSCC{2}),3)

L={}; L{2}=[1]; L{3}=[1]; L{4}=[1];
GSCC = tarjan(L);
assert(length(GSCC),4)



for iter=1:100  % random graph testing ....

  
  % testing undirected graphs
  adj = [0 1; 0 0];  % initialize so that the while loop does not break
  while not(isConnected(adj)); adj = randomGraph(randi(50)+1,rand); end

  % test that the GSCC contains one component, namely the entire graph
  L=adj2adjL(adj);
  GSCC = tarjan(L);
  assert(length(GSCC),1)
  assert(GSCC{1},[1:length(adj)])
  
  % testing directed graph
  adj=randomDirectedGraph(randi(50)+1,rand*0.1);
  L=adj2adjL(adj);
  GSCC = tarjan(L);
  
  
  if isConnected(adj) && isConnected(transpose(adj)) && length(adj)>0
    
    % there should be one component containing all nodes
    assert(length(GSCC),1)
    assert(GSCC{1},[1:length(adj)])
    
    
  else  % disconnected directed graph
    
    ll=[];
    for gg=1:length(GSCC); ll=[ll length(GSCC{gg})]; end;
    [ml,maxll]=max(ll);
    
    % the largest strongly connected component either is trivial (one node), or is connected
    assert(isConnected(adj(GSCC{maxll},GSCC{maxll})) | length(GSCC{maxll})==1)
    
    % adding any other node to the largest strongly connected component should not make the subgraph connected
    for ii=1:length(adj)
      if isempty(find(GSCC{maxll}==ii))
        
        tryGC = [GSCC{maxll}, ii];
        assert( not(isConnected(adj(tryGC,tryGC))) )
        
      end
      
    end
    
  end
  
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing subgraph.m ............................
fprintf('testing subgraph.m\n')
tic
assert(T{13}{2},subgraph(T{4}{2},[1,2,3]))
assert(T{13}{2},subgraph(T{4}{2},[4,5,6]))
assert(T{2}{2},subgraph(T{4}{2},[4,5]))
assert(T{2}{2},subgraph(T{4}{2},[1,2]))
assert(T{2}{2},subgraph(T{4}{2},[3,4]))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% ................................................
% ........ centrality measures ...................
% ................................................


% Testing sortNodesBySumNeighborDegrees.m ........
printf('testing sortNodesBySumNeighborDegrees.m\n')
tic
assert(sortNodesBySumNeighborDegrees(T{1}{2}),[1, 2]')
assert(sortNodesBySumNeighborDegrees(T{2}{2}),[2, 1]')
assert(sortNodesBySumNeighborDegrees(T{4}{2}),[4,3,6,5,2,1]')  
assert(sortNodesBySumNeighborDegrees(edgeL2adj(T{10}{2})),[1, 3, 2]')
assert(sortNodesBySumNeighborDegrees(T{13}{2}),[3, 2, 1]')
assert(sortNodesBySumNeighborDegrees(adjL2adj(T{17}{2})),[3, 2, 1]')
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing sortNodesByMaxNeighborDegree.m .........
printf('testing sortNodesByMaxNeighborDegree.m\n')
tic
assert(sortNodesByMaxNeighborDegree(T{2}{2}),[2, 1]')
assert(sortNodesByMaxNeighborDegree(T{4}{2}),[4,3,6,5,2,1]')  
assert(sortNodesByMaxNeighborDegree(edgeL2adj(T{10}{2})),[1, 3, 2]')
assert(sortNodesByMaxNeighborDegree(T{13}{2}),[3, 2, 1]')
assert(sortNodesByMaxNeighborDegree(adjL2adj(T{17}{2})),[3, 2, 1]')
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing transitivity.m .........................
printf('testing transitivity.m\n')
tic
assert(clustCoeff(T{13}{2}),1)
assert(clustCoeff(edgeL2adj(T{10}{2})),0)
assert(clustCoeff(edgeL2adj(canonicalNets(randi(10)+5,'tree',2))),0)
C = transitivity(T{4}{2});
assert(C,0.6)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% testing weightedClustCoeff.m ...................
printf('testing weightedClustCoeff.m\n')
tic
randint = randi(20);
assert(length(weightedClustCoeff(randomGraph(randint+5,rand))),randint+5)

adj = [0 2 1; 2 0 0; 1 0 0];
wC = weightedClustCoeff(adj);
assert(wC, [0 0 0]')

adj = [0 2 1; 2 0 1; 1 1 0];
wC = weightedClustCoeff(adj);
assert(wC, [1 1 1]')
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing sMetric.m ..............................
printf('testing sMetric.m\n')
assert(sMetric(T{13}{2}),2*12)
assert(sMetric(T{4}{2}),2*41)
assert(sMetric(edgeL2adj(T{11}{2})),4)
assert(sMetric(T{1}{2}),1)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% ................................................
% .............. distances .......................
% ................................................


% Testing simpleDijkstra.m .......................
printf('testing simpleDijkstra.m\n')
tic
assert(simpleDijkstra(T{4}{2},1),[0, 1, 1, 2, 3, 3])
assert(simpleDijkstra(T{4}{2},3),[1, 1, 0, 1, 2, 2])

mat = [0 3.5 0 1; 3.5 0 1 0; 0 1 0 1.4; 1 0 1.4 0];
assert(simpleDijkstra(mat,1),[0, 3.4, 2.4, 1])

assert(simpleDijkstra(edgeL2adj(T{11}{2}),1),[0, 1, 1])
assert(simpleDijkstra(edgeL2adj(T{11}{2}),2),[inf, 0, inf])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing dijkstra.m .............................
printf('testing dijkstra.m\n')
tic
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing shortestPathDP.m .......................
printf('testing shortestPathDP.m\n')
tic
bowtie = T{4}{2}
[Jb,rb,J,r]=shortestPathDP(T{4}{2},1,3,size(bowtie,1));
assert(Jb,1)
assert(rb,[1,3])

[Jb,rb,J,r]=shortestPathDP(T{4}{2},1,4,size(bowtie,1));
assert(Jb,2)
assert(rb,[1,3,4])

[Jb,rb,J,r]=shortestPathDP(T{4}{2},1,5,size(bowtie,1));
assert(Jb,3)
assert(rb,[1,3,4,5])

[Jb,rb,J,r]=shortestPathDP(edgeL2adj(T{11}{2}),1,2,3);
assert(Jb,1)
assert(rb,[1,2])

[Jb,rb,J,r]=shortestPathDP(edgeL2adj(T{11}{2}),2,3,3);
assert(Jb,inf)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing avePathLength.m ........................
printf('testing avePathLength.m\n')
tic
assert(avePathLength(T{4}{2}),(0+1+1+2+3+3 +0+1+2+3+3+ 0+1+2+2 +0+1+1 +0+1 +0)/15)
assert(avePathLength(T{13}{2}),1)
adj = edgeL2adj(canonicalNets(6,'line'));
assert(avePathLength(adj),(0+1+2+3+4+5 +0+1+2+3+4 +0+1+2+3 +0+1+2+ 0+1 +0)/15)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing smoothDiameter.m .......................
printf('testing smoothDiameter.m\n')
tic
adj = [0 1; 0 0];
while not(isConnected(adj)); adj = randomGraph(randi(10)+10,rand); end
assert(diameter(adj),smoothDiameter(adj,1))  % should be the same when the fraction is 1

assert( smoothDiameter(T{13}{2},1), 1 )
assert( smoothDiameter(T{13}{2},0.9999), 0 )
assert( smoothDiameter(T{13}{2},0.5), 0 )
assert( smoothDiameter(T{13}{2},0.1), 0 )
assert( smoothDiameter(T{13}{2},0), 0 )

assert( smoothDiameter(T{4}{2},1), 3 )
assert( smoothDiameter(T{4}{2}, 11/15), 2 ) % 7 node pairs at
                                            % diameter-2
assert( smoothDiameter(T{4}{2}, (7/15+11/15)/2), 1.5 )   % half
                                                        % between 1
                                                        % and 2
assert( smoothDiameter(T{4}{2}, 7/15), 1 )  % 7 node pairs at
                                            % diameter-1
assert( smoothDiameter(T{4}{2}, 6/15), 0 )
assert( smoothDiameter(T{4}{2}, 5/15), 0 )
assert( smoothDiameter(T{4}{2}, 0), 0 )
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing vertexEccentricity.m ...................
printf('testing vertexEccentricity.m\n')
tic
assert(vertexEccentricity(T{4}{2}),[3,3,2,2,3,3])
assert(vertexEccentricity(T{13}{2}),[1,1,1])
assert(vertexEccentricity(T{1}{2}),[1,inf])
assert(vertexEccentricity(edgeL2adj(T{11}{2})), [1, inf, inf])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% ................................................
% ....... simple motifs ..........................
% ................................................


% ................................................
% ......... linear algebra routines ..............
% ................................................


% Testing signlessLaplacian.m ......................
printf('testing signlessLaplacian.m\n')
tic
assert(signlessLaplacian(T{4}{2}),[2 1 1 0 0 0; 1 2 1 0 0 0; 1 1 3 1 0 0; 0 0 1 3 1 1; 0 0 0 1 2 1; 0 0 0 1 1 2])
assert(signlessLaplacian(T{13}{2}),[2 1 1; 1 2 1; 1 1 2])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% ................................................
% ......... modularity functions .................
% ................................................


% Testing simpleSpectralPartitioning.m ...........
printf('testing simpleSpectralPartitioning.m\n')
tic
for xx=1:50  % do the randomized test 50 times
  n = randi(99)+11;   % number of nodes
  adj = randomModularGraph(n,4,0.1,3);  % random graph with n nodes
  num_groups = randi(10)+1;  % number of groups to split the nodes in
  groups = [];
  for x=1:length(num_groups)-1; groups = [groups ceil(rand*n/num_groups)+1]; end
  groups = [groups n-sum(groups)];

  modules = simpleSpectralPartitioning(adj,groups);
  allnodes = [];
  for m=1:length(modules); 
      assert(length(modules{m}),groups(m)); 
      allnodes = [allnodes modules{m}];
  end
  assert( sort(allnodes), 1:n )
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................