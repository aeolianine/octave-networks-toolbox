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


% Testing rewire.m ...............................
printf('testing rewire.m\n')
tic
for x=1:100
  
  el = adj2edgeL(randomGraph(randi(10)+10,0.4));
  deg = degrees(edgeL2adj(el));
  rew = randi(5);
  eln = rewire(el,rew);
  degn = degrees(edgeL2adj(eln));
  
  assert(deg,degn)
  eq = eln(:,1:2) == el(:,1:2);
  preservedEdges = sum(sum(transpose(eq))==2);
  assert( preservedEdges >= size(eln)(1) - 4*rew )

end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing rewireThisEdge.m .......................
printf('testing rewireThisEdge.m\n')
tic
for x=1:100
  
  adj = [0 1; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(10)+10,0.3); end
  el = adj2edgeL(adj);
  deg = degrees(edgeL2adj(el));
  
  edgeind = randi([1,length(el)]);
  eln = rewireThisEdge(el,el(edgeind,1),el(edgeind,2));
  if isempty(eln); continue; end  % could not rewire, eln=[]
  
  adjn = edgeL2adj(eln);
  degn = degrees(adjn);
  
  assert(deg,degn)
  assert(isSimple(adjn),true)

  eq = eln(:,1:2) == el(:,1:2);
  preservedEdges = sum(sum(transpose(eq))==2);
  assert( preservedEdges == size(eln)(1) - 4 )

end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing rewireAssort.m .........................
printf('testing rewireAssort.m\n')
tic
for x=1:100
  adj = [0 0; 0 0];
  
  while not(isConnected(adj)); adj = randomGraph(randi(10)+10,0.4); end
  el = adj2edgeL(adj);
  eln = rewireAssort(el,randi(5));

  adjn = edgeL2adj(eln);
  assert(degrees(adj),degrees(adjn))

  assert(pearson(edgeL2adj(eln))>=(pearson(edgeL2adj(el))-10^(-7)) | (pearson(edgeL2adj(eln))-10^(-7))>=pearson(edgeL2adj(el)))
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing rewireDisassort.m ......................
printf('testing rewireDisassort.m\n')
tic
for x=1:100
  adj = [0 0; 0 0];
  
  while not(isConnected(adj)); adj = randomGraph(randi(10)+10,0.4); end
  el = adj2edgeL(adj);
  eln = rewireDisassort(el,randi(5));

  adjn = edgeL2adj(eln);
  assert(degrees(adj),degrees(adjn))

  assert(pearson(edgeL2adj(eln))<=(pearson(edgeL2adj(el))-10^(-7)) | (pearson(edgeL2adj(eln))-10^(-7))<=pearson(edgeL2adj(el)))
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
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


% Testing richClubMetric.m .......................
printf('testing richClubMetric.m\n')
assert(richClubMetric(randomGraph(randi(5)+5,rand),12),0)
assert(richClubMetric(T{4}{2},2),linkDensity(T{4}{2}))
assert(richClubMetric(T{4}{2},3),1)
assert(richClubMetric(T{4}{2},4),0)

mat = [0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0];
assert(richClubMetric(mat,2),1)
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
% ......... graph models .........................
% ................................................


% Testing randomGraph.m ..........................
printf('testing randomGraph.m\n');
tic
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing randomGraphFromDegreeSequence.m ........
printf('testing randomGraphFromDegreeSequence.m\n')
tic
for x=1:40
  
  adj = [0 1; 0 0];
  N = randi(50)+10;
  while not(isConnected(adj)); adj = randomGraph(N,log(N)/N); end
  
  adjr = randomGraphFromDegreeSequence(degrees(adj));
  
  assert(isSimple(adjr),true)
  assert(degrees(adj),degrees(adjr))
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing randomGraphDegreeDist.m ................
printf('testing randomGraphDegreeDist.m\n')
tic
pkg load statistics
N = randi(80)+30;
adj = randomGraphDegreeDist(N,'uniform');
assert(numNodes(adj),N)
assert(isSimple(adj),true)

[xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
plot(xpdf, ypdf, color = 'k')
text(5,0.01,strcat('constructing graphs with different degree distributions, N=  ', num2str(N)))
hold off; hold on;

adj = randomGraphDegreeDist(N,'normal');
assert(numNodes(adj),N)
assert(isSimple(adj),true)

[xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
plot(xpdf, ypdf, color = 'b')
hold off; hold on;
    
adj = randomGraphDegreeDist(N,'binomial');
assert(numNodes(adj),N)
assert(isSimple(adj),true)

[xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
plot(xpdf, ypdf, color = 'y')
hold off; hold on;


adj = randomGraphDegreeDist(N,'exponential');
assert(numNodes(adj),N)
assert(isSimple(adj),true)

[xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
plot(xpdf, ypdf, color = 'r')

legend('uniform','normal','binomial','exponential')

hold off;

adj = randomGraphDegreeDist(6,'custom',[1/5 1/5 1/5 1/5 1/5]);
assert(numNodes(adj),6)
assert(isSimple(adj),true)

adj = randomGraphDegreeDist(N,'custom');
assert(isempty(adj),true)

adj = randomGraphDegreeDist(N,'anything here');
assert(isempty(adj),true)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing randomModularGraph.m ...................
printf('testing randomModularGraph.m\n');
tic
for x=1:10
  N = randi(50)+10;
  c = randi(5)+1;
  [adj, modules] = randomModularGraph(N,c,0.2,4);
  assert(numNodes(adj),N)
  assert(length(modules),c)
end

% ..... testing with fixed labels ................
[adj, modules] = randomModularGraph(5,2,0.2,0.9,[1,1,2,3,4]);
mods = {}; mods{1} = [1,2]; mods{2} = [3]; mods{3} = [4]; mods{4} = 5;
assert( modules, mods )
[adj, modules] = randomModularGraph(6,0,0.2,0.8,[1,1,1,2,2,2]);
mods = {}; mods{1} = [1,2,3]; mods{2} = [4,5,6];
assert( modules, mods )
[adj, modules] = randomModularGraph(4,0,0.2,0.5,[1,2,3,4]);
mods = {}; mods{1} = 1; mods{2} = 2; mods{3} = 3; mods{4} = 4;
assert( modules, mods )
[adj, modules] = randomModularGraph(4,0,0.2,0.5,[1,1,1,1]);
mods = {}; mods{1} = [1,2,3,4];
assert( modules, mods )
[adj, modules] = randomModularGraph(75,4,0.2,0.5,[ones(1,20),2*ones(1,20),3*ones(1,20),4*ones(1,15)]);
mods = {}; mods{1} = [1:20]; mods{2} = [21:40]; mods{3} = [41:60]; mods{4} = [61:75];
assert( modules, mods )
% ................................................

% ... testing the in/out-degree ratio distribution ......

ratio = [];

for x=1:100
  
  mods = 4;
  N = 100;
  dens = log(N)/N;  % threshold of connectivity
  
  adj = [];
  while not(isConnected(adj)); [adj,modules] = randomModularGraph(N,mods,dens,10); end
  
  for M=1:mods

    kin = [];
    kout = [];
    
    for node=1:length(modules{M})
      i = modules{M}(node);

      ss = 0;
      for c=1:length(modules)
        if c == M; continue; end
        ss += sum(adj(i,modules{c}));
      end
  
      if ss==0; ss=0.01; end
      
      kin = [kin sum(adj(i,modules{M}))];
      kout = [kout ss];

      
    end
    
    kin = sum(kin)/length(kin);
    kout = sum(kout)/length(kout);
    
    ratio = [ratio kin/kout];
  end
 
  assert(length(modules),mods)
  assert(length(adj),N)
end

hist(ratio,20)
title('supposed to be centered around 10')
hold off;
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% .... end of test of randomModularGraph.m .......


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