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


% Testing numNodes.m ............................
printf('testing numNodes.m\n')
tic
randint = randi(101);
assert(numNodes(randomGraph(randint)),randint)

for i=1:length(T)
    if strcmp(T{i}{3},'adjacency')
        assert( numNodes(T{i}{2}), T{i}{6} )
    end
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing numEdges.m ...........................
printf('testing numEdges.m\n')
tic
for i=1:length(T)
    if strcmp(T{i}{3},'adjacency')
        assert( numEdges(T{i}{2}), T{i}{7} )
    end
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing linkDensity.m .........................
printf('testing linkDensity.m\n')
tic
randint = randi(101)+1;
assert(linkDensity(edgeL2adj(canonicalNets(randint,'tree',2))),2/randint)

for i=1:length(T)
    if strcmp(T{i}{3},'adjacency')
        coeff = 2;
        if isDirected(T{i}{2}); coeff = 1; end
        assert( linkDensity(T{i}{2}), coeff*T{i}{7}/(T{i}{6}*(T{i}{6}-1)) )
    end
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing selfLoops.m ...........................
printf('testing selfLoops.m\n')
tic
assert( selfLoops( edgeL2adj( T{8}{2} ) ), 1 )
assert( selfLoops( T{14}{2} ), 2 )
assert(selfLoops(T{4}{2}),0)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing multiEdges.m ..........................
printf('testing multiEdges.m\n')
tic
assert(multiEdges(T{3}{2}),1)
assert(multiEdges([0 2 1; 2 0 1; 1 1 0]),1)  % triangle with one double edge
assert(multiEdges([0 0 1; 2 0 0; 0 1 0]),1)  % directed triangle with 1 double edge
assert(multiEdges(randomGraph(randi(15))),0)
assert(multiEdges([0 0 1; 2 0 0; 0 2 0]),2)  % directed triangle with 2 double edges
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing numConnComp.m .........................
printf('testing numConnComp.m\n')
tic
assert(numConnComp(T{5}{2}),2)

randint = randi(51);
Adj=zeros(randint*30);
for x=1:randint
  adj=randomGraph(30,0.5);
  Adj(30*(x-1)+1:30*x,30*(x-1)+1:30*x)=adj;
end
assert(numConnComp(Adj),randint)
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

% Testing leafNodes.m ...........................
fprintf('testing leafNodes.m\n')
tic
assert(leafNodes(edgeL2adj(T{10}{2})),[2,3])
assert(leafNodes(edgeL2adj(T{11}{2})),[2,3])
assert(length(leafNodes(T{13}{2})),0)
assert(leafNodes(T{2}{2}),[1,2])
assert(leafNodes(T{1}{2}),[2])
assert(length(leafNodes(T{4}{2})),0)
assert(leafNodes(edgeL2adj(T{19}{2})),[2,3,4,5])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................

% Testing leafEdges.m ...........................
fprintf('testing leafEdges.m\n')
tic
assert(leafEdges(edgeL2adj(T{10}{2})),[1,2;1,3])
assert(leafEdges(edgeL2adj(T{10}{2})),[1,2;1,3])
assert(length(leafEdges(T{13}{2})),0)
assert(length(leafEdges(edgeL2adj([2,1,1;3,1,1]))),0)
assert(length(leafEdges(T{4}{2})),0)
assert(leafEdges(edgeL2adj(T{19}{2})),[1, 2; 1, 3; 1, 4; 1, 5])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................

% Testing minSpanTree.m .........................
printf('testing minSpanTree.m\n')
tic
for x=1:100
  adj = [0 1; 0 0];                % initialize
  while not(isConnected(adj)); adj = randomGraph(randi(30)+5,rand); end
  
  tr = minSpanTree(adj);
  assert(isTree(tr),true)
  assert(length(tr),length(adj));  % tree should have the same
                                   % number of nodes as adj
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................

% ................................................
% ........ diagnostic functions ..................
% ................................................


% Testing isSimple.m .............................
fprintf('testing isSimple.m\n')
tic
assert(isSimple(T{1}{2}),false)
assert(isSimple(T{2}{2}),true)
assert(isSimple(T{3}{2}),false)
assert(isSimple(randomGraph(randi(5)+20,rand)),true)  % simple graph
assert(isSimple(edgeL2adj([1,2,2])),false)      % multi-edge
assert(isSimple( [1 0 0; 0 0 1; 0 1 0]),false)  % matrix with loops
assert(isSimple([0 1 1; 1 0 0; 0 1 0]),false)   % directed matrix
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing isDirected.m ...........................
fprintf('testing isDirected.m\n')
tic
assert(isDirected(randomDirectedGraph(randi(5)+20,rand)),true)  
assert(isDirected(randomGraph(randi(5)+20,rand)),false)
assert(isDirected(T{1}{2}), true)
assert(isDirected(T{2}{2}), false)
assert(isDirected(T{3}{2}), false)
assert(isDirected(T{16}{2}), true)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing isSymmetric.m ..........................
fprintf('testing isSymmetric.m\n')
tic
for i=1:100
  assert(isSymmetric(randomGraph(randi(5)+20,rand)),true)

  adj = randomDirectedGraph(randi(5)+20,rand);
  assert(not(isSymmetric(adj)) | adj==zeros(size(adj)) | adj==ones(size(adj)))
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing isConnected.m ..........................
fprintf('testing isConnected.m\n')
tic
assert(isConnected(T{1}{2}),false)
assert(isConnected(T{2}{2}),true)
assert(isConnected(T{4}{2}),true)
assert(isConnected(T{5}{2}),false)
assert(isConnected(transpose(T{5}{2})),false)
assert(isConnected(T{13}{2}),true)
assert(isConnected(T{14}{2}),true)

for x=1:100
  % test the connected case
  adj = [0 1; 1 0];                % initialize
  N = randi(100)+5;
  while isConnected(adj) || sum(adj)==0; adj = randomGraph(N,log(N)/N); end
  assert(isConnected(adj),false)
  assert(abs(algebraicConnectivity(adj)) <10^(-10),true)
  assert(transpose(isConnected(adj)),false)
  assert(length(findConnComp(adj))>1, true)
  
  % test the not-connected case
  adj = [0 1; 0 0];                % initialize
  while not(isConnected(adj)); adj = randomGraph(N,log(N)/N); end
  assert(isConnected(adj),true)
  assert(algebraicConnectivity(adj)>0,true)
  assert(transpose(isConnected(adj)),true)
  assert(length(findConnComp(adj))==1, true)
  
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing isWeighted.m ...........................
fprintf('testing isWeighted.m\n')
tic
assert(isWeighted(adj2edgeL(T{2}{2})),false)
assert(isWeighted(adj2edgeL(T{3}{2})),true)

assert(isWeighted(adj2edgeL(randomGraph(randi(5)+20,rand+0.1))),false)
  
assert(isWeighted(adj2edgeL(randomDirectedGraph(randi(5)+20,rand+0.1))),false)
  
assert(isWeighted([1,2,0.5; 1,3,1.5; 1,4,1]),true)
assert(isWeighted([1,2,0.5; 1,3,1; 1,4,1]),true)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing isRegular.m ............................
fprintf('testing isRegular.m\n')
tic
adj = edgeL2adj(canonicalNets(20,'cycle'));
assert(isRegular(adj),true)

adj = edgeL2adj(canonicalNets(20,'tree',3));
assert(isRegular(adj),false)

assert(isRegular([0 1; 1 0]),true)
assert(isRegular([0 0; 1 0]),false)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing isComplete.m ...........................
printf('testing isComplete.m\n')
tic
assert(isComplete([0 1; 1 0]),true)
assert(isComplete(T{2}{2}),true)
assert(isComplete(T{3}{2}),true)
assert(isComplete(T{4}{2}),false)

randint = randi(10)+10;
adj = ones(randint)-eye(randint);
assert(isComplete(adj),true)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing isEulerian.m ...........................
printf('testing isEulerian.m\n')
tic
adj = edgeL2adj(canonicalNets(10,'cycle'));
assert(isEulerian(adj),true)

adj = edgeL2adj(canonicalNets(10,'tree',3));
assert(isEulerian(adj),false)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing isTree.m ...............................
printf('testing isTree.m\n')
tic
assert(isTree(T{1}{2}), false)
assert(isTree(T{2}{2}), true)
assert(isTree(T{3}{2}), false)
assert(isTree(T{4}{2}), false)
assert(isTree(T{5}{2}), false)
assert(isTree(edgeL2adj(T{10}{2})), true)
assert(isTree(edgeL2adj(T{11}{2})), false)
assert(isTree(T{13}{2}), false)
assert(isTree(T{14}{2}), false)
assert(isTree(T{16}{2}), false)
assert(isTree(T{18}{2}), false)
assert(isTree(edgeL2adj(T{19}{2})), true)

adj = edgeL2adj(canonicalNets(randi(10)+10,'hexlattice'));
assert(isTree(adj),false)
adj = edgeL2adj(canonicalNets(randi(10)+10,'trilattice'));
assert(isTree(adj),false)
adj = edgeL2adj(canonicalNets(randi(10)+10,'hierarchy', b=3));
assert(isTree(adj),false)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing isGraphic.m ............................
printf('testing isGraphic.m\n')
tic
for i=1:20
  adj = giantComponent(randomGraph(randi(20)+1,0.5));
  [deg,~,~] = degrees(adj);
  assert(isGraphic(deg) | adj==0)
end

assert(isGraphic([0 1]), false)
assert(isGraphic([2 1]), false)
assert(isGraphic([1 1 4]), false)
assert(isGraphic([1 4 4 100]), false)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing isBipartite.m ..........................
printf('testing isBipartite.m\n')
tic

% ................................................

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

% Testing nodeBetweenness.m ......................
printf('testing nodeBetweenness.m and nodeBetweennessFaster\n')
tic
assert(nodeBetweenness([0 1; 1 0]),[0 0])
assert(nodeBetweennessFaster([0 1; 1 0]),[0 0])

assert(nodeBetweenness([1 1; 0 0]),[0 0])
assert(nodeBetweennessFaster([1 1; 0 0]),[0 0])

assert(nodeBetweenness([0 1 1; 1 0 0; 1 0 0]),[1/3 0 0])
assert(nodeBetweennessFaster([0 1 1; 1 0 0; 1 0 0]),[1/3 0 0])

assert(nodeBetweenness(T{4}{2}),[0 0 0.4 0.4 0 0])
assert(nodeBetweennessFaster(T{4}{2}),[0 0 0.4 0.4 0 0])

x = edgeL2adj(canonicalNets(2*randi(10)+2,'cycle'));
bw = nodeBetweenness(x);
assert(bw(1)*ones(1,length(bw)),bw)  % the betweennesses should be all the same
bw = nodeBetweennessFaster(x);
assert(bw(1)*ones(1,length(bw)),bw)  % the betweennesses should be all the same

L={}; L{1}=[2]; L{2}=[1,3,4]; L{3}=[2,5]; L{4}=[2,5]; L{5}=[3,4,6]; L{6}=[5];
adj = adjL2adj(L);
bw = nodeBetweenness(adj);
bwF = nodeBetweennessFaster(adj);
assert(bw,bwF)

adj = [];
while not(isConnected(adj)); adj = randomGraph(20,log(20)/20); end
bw = nodeBetweenness(adj);
bwF = nodeBetweennessFaster(adj);
assert(norm(bw-bwF)<10^(-10))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing edgeBetweenness.m ......................
printf('testing edgeBetweenness.m\n')
tic
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


% Testing pearson.m ..............................
printf('testing pearson.m\n')
tic
assert(pearson( edgeL2adj(T{10}{2}) ),-1)
assert(pearson( edgeL2adj(T{19}{2}) ),-1)
assert( pearson( edgeL2adj( canonicalNets(randi(5)+5,'star') ) ) ,-1 )

% test via pearsonW.m (Whitney routine)
for i=1:50
  if isComplete(adj); continue; end
  adj = randomGraph(randi(20)+3,rand);
  assert( abs(pearson(adj)-pearsonW(adj))<10^(-6)  )
end
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


% Testing kneighbors.m ...........................
printf('testing kneighbors.m\n')
tic
assert(kneighbors(T{4}{2},1,3),[1 2 3 4 5 6])
assert(kneighbors(T{4}{2},3,1),[1 2 4])
assert(kneighbors(T{13}{2},2,1),[1,3])
assert(kneighbors(T{13}{2},1,2),[1,2,3])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing kminNeighbors.m ........................
fprintf('testing kminNeighbors.m\n')
tic
assert(kminNeighbors(T{4}{2},1,3),[5, 6])
assert(kminNeighbors(T{4}{2},3,1),[1, 2, 4])
assert(kminNeighbors(T{4}{2},3,2),[5, 6])
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


% Testing numConnTriples.m .......................
printf('testing numConnTriples.m\n')
tic
assert(numConnTriples(T{4}{2}),6)
assert(numConnTriples(T{13}{2}),1)
assert(numConnTriples(edgeL2adj(T{10}{2})),1)
assert(numConnTriples(T{2}{2}),0)
assert(numConnTriples(T{18}{2}),4)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing numCycles.m .............................
printf('testing numCycles.m\n')
tic
assert(numCycles(T{13}{2}),1)
assert(numCycles(T{4}{2}),2)
assert(numCycles(edgeL2adj(T{10}{2})),0)
assert(numCycles(T{18}{2}),1)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing loops3rev2.m ...............................
printf('testing loops3rev2.m\n')
tic
fourCycle = T{18}{2};
assert(loops3rev2(T{4}{2}),{'1-2-3','4-5-6'})
assert(loops3rev2(fourCycle),{})
assert(loops3rev2(T{13}{2}),{'1-2-3'})
assert(loops3rev2(edgeL2adj(canonicalNets(randi(10)+3,'btree'))),{})
assert(loops3rev2(edgeL2adj(canonicalNets(4,'trilattice'))),{'1-2-4','1-3-4'})
assert(loops3rev2(T{16}{2}),{'1-2-3'})


printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing numStarMotifs.m ........................
printf('testing numStarMotifs.m\n')
tic
assert(numStarMotifs(T{9}{2},3),4+6)
assert(numStarMotifs(T{9}{2},4),2)
assert(numStarMotifs(T{9}{2},5),0)

assert(numStarMotifs(adj2adjL(T{13}{2}),3),3)
assert(numStarMotifs(adj2adjL(T{13}{2}),2),6)

assert(numStarMotifs(T{9}{2},1),6)   % trivial case
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% ................................................
% ......... linear algebra routines ..............
% ................................................

% Testing laplacianMatrix.m ......................
printf('testing laplacianMatrix.m\n')
tic
assert(laplacianMatrix(T{4}{2}),[2 -1 -1 0 0 0; -1 2 -1 0 0 0; -1 -1 3 -1 0 0; 0 0 -1 3 -1 -1; 0 0 0 -1 2 -1; 0 0 0 -1 -1 2])
assert(laplacianMatrix(T{13}{2}),[2 -1 -1; -1 2 -1; -1 -1 2])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
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


% test kregular.m ................................
printf('testing kregular.m\n');
tic
for x=1:30
  
  n = randi(20)+5;   % random integer between 6 and 25
  k = randi(n-2)+1;  % randon integer between 2 and n-1
  if mod(k,2)==1 && mod(n,2)==1; continue; end  % no solution for this case
  el = kregular(n,k);
  adj = edgeL2adj(el);
  assert(degrees(adj),k*ones(1,length(adj)))

end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing randomGraph.m ..........................
printf('testing randomGraph.m\n');
tic
% testing the size of the graph
randint = randi(20)+3;
assert(size(randomGraph(randint),1),randint)
assert(size(randomGraph(randint),2),randint)

% testing the default probability of attachment
for x=1:50
  randint = randi(50)+50;
  adj = randomGraph(randint);
  assert(linkDensity(adj)>0.4);
  assert(linkDensity(adj)<0.6);
end

% testing a random probability of attachment
for x=1:50
  p = rand;
  randint = randi(50)+50;
  adj = randomGraph(randint,p);
  assert(linkDensity(adj)>p-0.05);
  assert(linkDensity(adj)<p+0.05);
end

% testing for the number of edges, E
for x=1:50
  randint = randi(50)+50;
  E = randi([1,randint-1]);
  adj = randomGraph(randint,[],E);
  assert(numEdges(adj),E);
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing randomDirectedGraph.m ..................
printf('testing randomDirectedGraph.m\n');
tic
for i=1:30
  p=rand;
  n = randi(40)+40;
  adj = randomDirectedGraph(n,p);
  assert(linkDensity(adj)>p-0.05)
  assert(linkDensity(adj)<p+0.05)
  assert(sum(adj)==0 || isDirected(adj))
  assert(size(adj),[n,n]);
end
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


% Testing PriceModel.m ...........................
printf('testing PriceModel.m\n')
tic
for x=1:20
  randint = randi(10)+10;
  adj = PriceModel(randint);
  assert(isDirected(adj),true)
  assert(numNodes(adj),randint)    
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing preferentialAttachment.m ...............
printf('testing preferentialAttachment.m\n')
tic
for x=1:10
  el = preferentialAttachment(randi(10)+10,1);
  adj = edgeL2adj(el);
  assert(isTree(adj),true)
  assert(isSimple(adj),true)
  
  randint = randi(30)+5;
  el = preferentialAttachment(randint,2);
  adj = edgeL2adj(el);
  assert(numEdges(adj),1+2*(length(adj)-2))
    
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing masterEquation.m .......................
printf('testing masterEquation.m\n')
tic
for x=1:30
  randint = randi(100)+5;
  adj = masterEquationGrowthModel(randint,1,0);
  assert(isTree(adj),true)
  
  adj = masterEquationGrowthModel(randint,2,0);
  assert(isTree(adj),false)
  
  adj = masterEquationGrowthModel(randint,2,2);
  assert(isSimple(adj),true)
  
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing newmanGastner.m ........................
printf('testing newmanGastner.m\n')
tic
for x=1:10
  N = randi(100)+10;
  el = newmanGastner(N,rand,[],'off');  % no plot
  adj = symmetrize(edgeL2adj(el));
  assert(numNodes(adj),N);
  assert(isSimple(adj),true)
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing nestedHierarchiesModel.m ...............
printf('testing nestedHierarchiesModel.m\n')
tic
el = nestedHierarchiesModel(640,3,[10, 20, 40],10);
adj = edgeL2adj(el);
assert(isSimple(adj));
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

% Testing newmanGirvan.m .........................
printf('testing newmanGirvan.m\n')
tic
[modules, moduleHist, Q] = newmanGirvan(T{4}{2},2);
assert(modules{1}==[1,2,3])
assert(modules{2}==[4,5,6])
assert(moduleHist{1}==[1,2,3,4,5,6])
assert(moduleHist{2}==[1,2,3])
assert(abs(Q-0.20408)<10^(-5))

[modules, moduleHist, Q] = newmanGirvan(edgeL2adj(T{19}{2}),2);
assert(modules{1}==[1,3,4,5])
assert(modules{2}==[2])
assert(moduleHist{1}==[1,2,3,4,5])
assert(moduleHist{2}==[1,3,4,5])
assert(abs(Q+0.31250)<10^(-5))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing newmanEigenvectorMethod.m ..............
printf('testing newmanEigenvectorMethod.m\n')
tic
modules = newmanEigenvectorMethod(T{4}{2});
assert(length(modules),2)
assert(modules{1},[4,5,6])
assert(modules{2},[1,2,3])

for x=1:50
  adj = randomGraph(randi(10)+5,1);
  Adj = zeros(4*length(adj));
  Adj(1:length(adj),1:length(adj))=adj;
  Adj(length(adj)+1:2*length(adj),length(adj)+1:2*length(adj))=adj;
  Adj(2*length(adj)+1:3*length(adj),2*length(adj)+1:3*length(adj))=adj;
  Adj(3*length(adj)+1:4*length(adj),3*length(adj)+1:4*length(adj))=adj;

  Adj(5,length(adj)+5)=1; Adj(length(adj)+5,5)=1; 
  Adj(length(adj)+6,2*length(adj)+6)=1; Adj(2*length(adj)+6,length(adj)+6)=1; 
  Adj(2*length(adj)+7,3*length(adj)+7)=1; Adj(3*length(adj)+7,2*length(adj)+7)=1; 
  Adj(3*length(adj)+1,1)=1; Adj(1,3*length(adj)+1)=1; 

  modules = newmanEigenvectorMethod(Adj);
  assert(length(modules),4)

  prescribed = randi(6)+2;
  
  n = randi(50)+50;
  adj = [];
  while not(isConnected(adj)); adj = randomModularGraph(n,prescribed,0.9*log(n)/n,2+0.3*rand); end
  modules = newmanEigenvectorMethod(adj);
  
  sumnodes = 0;
  for m=1:length(modules); sumnodes = sumnodes + length(modules{m}); end
  assert(sumnodes,n)
  
  for m1=1:length(modules)
    for m2=m1+1:length(modules)
      assert(length(intersect(modules{m1},modules{m2})),0)     
    end
  end
   
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing newmanCommFast.m .......................
printf('testing newmanCommFast.m\n')
tic
[gH,Q]=newmanCommFast(T{4}{2});
close all;
assert(max(Q),Q(6-1));

[gH,Q]=newmanCommFast(randomModularGraph(100,4,0.1,5));
close all;
assert(length(gH),length(Q))
[~,ind]=max(Q);
assert(length(gH{ind}),4)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing modularityMetric.m .....................
printf('testing modularityMetric.m\n')
tic
for i=1:10
  
  adj = [0 1; 0 0];
  num_modules = randi([2,5]);
  while not(isConnected(adj)); adj = randomModularGraph(30,num_modules,0.1,3); end
 
  % compare to newmanCommFast
  [mH,Q1] = newmanCommFast(adj);
  close all;
  Q2 = [];
  for m=1:length(mH); Q2 = [Q2 modularityMetric(mH{m},adj)]; end
  
  assert(Q1,Q2)

  % compare to the newman-girvan routine
  [modules0,~,Q0] = newmanGirvan(adj,num_modules);
  assert(Q0,modularityMetric(modules0,adj))
 
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing louvainCommunityFinding.m ..............
printf('testing louvainCommunityFinding.m\n');
tic
extended_bowtie0 = [1 2; 2 3; 3 2; 3 4; 4 5; 5 6; 4 6; 6 7; 7 8; 7 9; 8 9];
extended_bowtie = [];
for row=1:size(extended_bowtie0,1)
  extended_bowtie = [extended_bowtie; extended_bowtie0(row,:) 1];
end
clear extended_bowtie0
extended_bowtie = symmetrizeEdgeL(extended_bowtie);
adj = edgeL2adj(extended_bowtie);

[modules,inmodule]=louvainCommunityFinding(adj);
assert(length(modules),3)
assert([inmodule{1},inmodule{2},inmodule{3}],[1,1,1]*inmodule{1})
assert([inmodule{4},inmodule{5},inmodule{6}],[1,1,1]*inmodule{4})
assert([inmodule{7},inmodule{8},inmodule{9}],[1,1,1]*inmodule{7})

[modules,inmodule]=louvainCommunityFinding(bowtie);
assert(length(modules),2)
assert([inmodule{1},inmodule{2},inmodule{3}],[1,1,1]*inmodule{1})
assert([inmodule{4},inmodule{5},inmodule{6}],[1,1,1]*inmodule{4})

% concatenate 4 complete graphs: 
adj = ones(10,10)-eye(10);
Adj = zeros(40,40);

Adj(1:10,1:10)=adj;
Adj(11:20,11:20)=adj;
Adj(21:30,21:30)=adj;
Adj(31:40,31:40)=adj;

Adj(10,11) = 1; Adj(11,10) = 1;
Adj(20,21) = 1; Adj(21,20) = 1;
Adj(30,31) = 1; Adj(31,30) = 1;

[modules,inmodule]=louvainCommunityFinding(Adj);
assert(length(modules),4)

% concatenate 4 dense graphs
adj = [0 1; 0 0];
while not(isConnected(adj)); adj = randomGraph(10,0.9); end
Adj = zeros(40,40);

Adj(1:10,1:10)=adj;
Adj(11:20,11:20)=adj;
Adj(21:30,21:30)=adj;
Adj(31:40,31:40)=adj;

Adj(10,11) = 1; Adj(11,10) = 1;
Adj(20,21) = 1; Adj(21,20) = 1;
Adj(30,31) = 1; Adj(31,30) = 1;

[modules,inmodule]=louvainCommunityFinding(Adj);
assert(length(modules),4)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% ................................................
% ......... simple matrix/graph viz ..............
% ................................................


% Testing pdfCdfRank.m ...........................
printf('testing pdfCdfRank.m\n');
tic
adj = randomGraph(randi(30)+30,0.2);
[xp,yp,xc,yc,lk,lx] = pdfCdfRank(degrees(adj),'off');
assert(length(xp),length(xc))
assert(length(xp),length(yp))
assert(length(yp),length(yc))
assert(length(lk),length(lx))
[xp,yp,xc,yc,lk,lx] = pdfCdfRank(degrees(adj),'off', bin=5);
assert(length(xp),length(xc))
assert(length(xp),length(yp))
assert(length(yp),length(yc))
assert(length(lk),length(lx))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................
