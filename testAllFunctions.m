% Test code for "Octave routines for Network Analysis"

clear all
close all


% Set of test graphs, in various formats =========
one_double_edge = [0 2; 2 0]; 
bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];   % 'adj'
disconnected_bowtie =[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];   % 'adj'
bowtie_edgeL = [1,2,1; 1,3,1; 2,3,1; 3,4,1; 4,5,1; 4,6,1; 5,6,1]; % 'edgelist'
bowtie_edgeL = sortrows(symmetrizeEdgeL(bowtie_edgeL));          % 'edgelist'
bowtie_edgeL_loop = [bowtie_edgeL; 4 4 1];                        % 'edgelist'
bowtie_adjL = {[2,3],[1,3],[1,2,4],[3,5,6],[4,6],[4,5]};          % 'adjlist'
undirected_cherry = [1,2,1; 2,1,1; 1,3,1; 3,1,1];                 %
                                                                  % 'edgelist' undirected
directed_cherry = [1,2,1; 1,3,1];                                 % 'edgelist'
undirected_triangle=[0 1 1; 1 0 1; 1 1 0];                        % 'adj'
undirected_triangle_selfloops = [1 1 1; 1 1 1; 1 1 0];            % 'adj'
undirected_triangle_incidence = [1 1 0; 1 0 1; 0 1 1];            % 'adj'
directed_triangle=[0 1 0; 0 0 1; 1 0 0];                          % 'adj'
square = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];                    % 'adj'
% ================================================


% Testing getNodes.m =============================
printf('testing getNodes.m\n')

assert(getNodes(bowtie,'adj'), [1:6])
N = randi(100);
assert(getNodes(random_directed_graph(N),'adj'),[1:N])
assert(getNodes(random_graph(10),'adj'),[1:10])
assert(getNodes(bowtie_adjL,'adjlist'),[1:6])
assert(getNodes(directed_cherry,'edgelist'),[1:3])
assert(getNodes(undirected_cherry,'edgelist'),[1:3])
assert(getNodes(undirected_triangle_incidence,'inc'),[1:3])
% ================================================

% Testing getEdges.m =============================
printf('testing getEdges.m\n')

assert(getEdges(bowtie,'adj'),bowtie_edgeL)
assert(getEdges(bowtie_adjL,'adjlist'),bowtie_edgeL)
assert(getEdges(directed_cherry,'edgelist'),directed_cherry)
assert(getEdges(undirected_cherry,'edgelist'),undirected_cherry)
assert(getEdges(undirected_triangle_incidence,'inc'),[1,2,1; 1,3,1; 2,1,1; 2,3,1; 3,1,1; 3,2,1])
% ================================================

% testing numNodes.m =============================
printf('testing numNodes.m\n')

randint = randi(101);
assert(numNodes(random_graph(randint)),randint)
assert(numEdges(edgeL2adj(directed_cherry)),2)
assert(numNodes(bowtie),6)
% ================================================

% testing numEdges.m =============================
printf('testing numEdges.m\n')

assert(numEdges(bowtie),7)
assert( numEdges(undirected_triangle_selfloops), 5 )
assert(numEdges(one_double_edge),2)
assert(numEdges(edgeL2adj(bowtie_edgeL_loop)),8)
% ================================================


% testing linkDensity.m =========================
printf('testing linkDensity.m\n')

randint = randi(101);
assert(linkDensity(edgeL2adj(canonical_nets(randint,'tree',2))),2/randint)
assert(linkDensity(bowtie),2.0*7/(6*5))
% ================================================


% testing selfLoops.m ============================
printf('testing selfLoops.m\n')

assert(selfLoops(undirected_triangle_selfloops),2)
assert(selfLoops(bowtie),0)
% ================================================


% testing multiEdges.m ===========================
printf('testing multiEdges.m\n')

assert(multiEdges(one_double_edge),2)
assert(multiEdges([0 2 1; 2 0 1; 1 1 0],2))  % triangle with one double edge
assert(multiEdges([0 0 1; 2 0 0; 0 1 0]),2)  % directed triangle with 1 double edge
assert(multiEdges(random_graph(randi(15))),0)
% ================================================


% testing averageDegree.m ========================
printf('testing averageDegree.m\n')

assert(averageDegree(square),2)
assert(averageDegree(bowtie),2+1.0/3)
% ================================================


% testing numConnComp.m ==========================
printf('testing numConnComp.m\n')
nc=numConnComp(disconnected_bowtie);
assert(numConnComp(disconnected_bowtie),2)

randint = randi(51);
Adj=zeros(randint*30);
for x=1:randint
  adj=random_graph(30,0.5);
  Adj(30*(x-1)+1:30*x,30*(x-1)+1:30*x)=adj;
end
assert(numConnComp(Adj),randint)
% ================================================


% testing findConnComp.m =========================
printf('testing findConnComp.m\n')

assert(findConnCompI(disconnected_bowtie,1),[1,2,3])
assert(findConnComp(disconnected_bowtie),{[1,2,3],[4,5,6]})

clear modules
modules{1}=[0];
randint = randi(21);
Adj = []; adj = [];

% make up a matrix (Adj) of randint disconnected components (adj)
for x=1:randint
  randsecint = randi(25)+5;
  lastnode = modules{length(modules)}(length(modules{length(modules)}));
  modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
  while isempty(adj) | not(isConnected(adj)) | not(length(adj)==randsecint); adj=random_graph(randsecint,0.5); end

  Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
end

modules=modules(2:length(modules));
assert(findConnComp(Adj),modules)
% ================================================


% testing giantComponent.m =======================
printf('testing giantComponent.m\n')

clear modules
modules{1}=[0];
randint = randi(20)+1;
Adj = []; adj = [];

% make up a matrix (Adj) of randint disconnected components (adj)
for x=1:randint
  randsecint = randi(25)+5;
  lastnode = modules{length(modules)}(length(modules{length(modules)}));
  modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
  while isempty(adj) | not(isConnected(adj)) | not(length(adj)==randsecint); adj=random_graph(randsecint,0.5); end
  Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
end
modules=modules(2:length(modules));
L = [];
for m=1:length(modules); L = [L, length(modules{m})]; end;
[maxL,maxind] = max(L);
assert(giantComponent(Adj), subgraph(Adj,modules{maxind}))
% ================================================



% ================================================
% Testing tarjan.m ===============================
printf('testing tarjan.m\n')

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


for iter=1:100  % completely random matrix testing ....

  
  % undirected graph testing ========================
  adj = [0 1; 0 0];  % initialize so that the while loop does not break
  while not(isConnected(adj)); adj = random_graph(randi(50)+1,rand); end

  L=adj2adjL(adj);
  GSCC = tarjan(L);
  assert(length(GSCC),1)
  assert(GSCC{1},[1:length(adj)])
  
  % directed graph testing ==========================
  adj=random_directed_graph(randi(50)+1,rand);
  L=adj2adjL(adj);
  GSCC = tarjan(L);
  
  
  if isConnected(adj) & isConnected(transpose(adj)) & length(adj)>0
    
    % there should be one component containing all nodes
    assert(length(GSCC),1)
    assert(GSCC{1},[1:length(adj)])
    
    
  else  % disconnected directed graph
    
    ll=[];
    for gg=1:length(GSCC); ll=[ll length(GSCC{gg})]; end;
    [ml,maxll]=max(ll);
    
    assert(isConnected(adj(GSCC{maxll},GSCC{maxll})) | length(GSCC{maxll})==1)
    
    for ii=1:length(adj)
      if isempty(find(GSCC{maxll}==ii))
        
        tryGC = [GSCC{maxll}, ii];
        assert(not(isConnected(adj(tryGC,tryGC))) | not(isConnected(transpose(adj(tryGC,tryGC)))))
        
      end
      
    end
    
  end
  
end
% ================================================
% ================================================


% testing graphComplement.m =====================
printf('testing graphComplement.m\n')

mat = [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 0 1 1; 1 1 0 1 0 0; 1 1 1 0 1 0; 1 1 1 0 0 1];
assert(graphComplement(bowtie),mat)
assert(graphComplement(undirected_triangle),eye(3))  
% ================================================


% Testing graphDual.m ============================
printf('testing graphDual.m\n')

gd=graphDual(adj2adjL(bowtie));
gdT={};
gdT{1}=[2,3]; gdT{2}=[1,3,4]; gdT{3}=[1,2,4]; gdT{4}=[2,3,5,6]; gdT{5}=[4,6,7]; gdT{6}=[4,5,7]; gdT{7}=[5,6];
assert(gd,gdT)

gd=graphDual(adj2adjL(undirected_triangle));
gdT={};
gdT{1}=[2,3]; gdT{2}=[1,3]; gdT{3}=[1,2];
assert(gd,gdT)

L={}; LT={}; L{1}=[2]; L{2}=[1]; LT{1}=[];
assert(LT,graphDual(L))
% ================================================

% testing subgraph.m =============================
printf('testing subgraph.m\n')
assert(undirected_triangle,subgraph(bowtie,[1,2,3]))
% ================================================

% testing leafNodes.m ===========================
printf('testing leafNodes.m\n')
assert(leafNodes(edgeL2adj(undirected_cherry)),[2,3])
assert(leafNodes(edgeL2adj(directed_cherry)),[2,3])
assert(length(leafNodes(undirected_triangle)),0)
% ================================================

% testing leafEdges.m ===========================
printf('testing leafEdges.m\n')
assert(leafEdges(edgeL2adj(undirected_cherry)),[1,2;1,3])
assert(leafEdges(edgeL2adj(directed_cherry)),[1,2;1,3])
assert(length(leafEdges(undirected_triangle)),0)
hut = [2,1,1;3,1,1];
assert(length(leafEdges(edgeL2adj(hut))),0)
% ================================================

% testing issimple.m =============================
printf('testing isSimple.m\n')

assert(isSimple(random_graph(randi(5)+20,rand)),true)  % simple graph
assert(isSimple(edgeL2adj([1,2,2])),false)      % multi-edge
assert(isSimple( [1 0 0; 0 0 1; 0 1 0]),false)  % matrix with loops
assert(isSimple([0 1 1; 1 0 0; 0 1 0]),false)   % directed matrix
% ================================================

% testing isDirected.m ===========================
printf('testing isDirected.m\n')
assert(isDirected(random_directed_graph(randi(5)+20,rand)),true)  
assert(isDirected(random_graph(randi(5)+20,rand)),false)
% ================================================


% testing isSymmetric.m ==========================
printf('testing isSymmetric.m\n')

for i=1:100
  assert(isSymmetric(random_graph(randi(5)+20,rand)),true)

  adj = random_directed_graph(randi(5)+20,rand);
  assert(not(isSymmetric(adj)) | adj==zeros(size(adj)) | adj==ones(size(adj)))
end
% ================================================

% testing isConnected.m ==========================
printf('testing isConnected.m\n')
assert(isConnected(bowtie),true)
assert(isConnected(disconnected_bowtie),false)
% ================================================

% testing isWeighted.m ===========================
printf('testing isWeighted.m\n')
assert(isWeighted([1,2,2]),true)

assert(isWeighted(adj2edgeL(random_graph(randi(5)+20,rand))),false)
  
assert(isWeighted(adj2edgeL(random_directed_graph(randi(5)+20,rand))),false)
  
assert(isWeighted([1,2,0.5; 1,3,1.5; 1,4,1]),true)
assert(isWeighted([1,2,0.5; 1,3,1; 1,4,1]),true)
% ================================================

% testing isRegular.m ============================
printf('testing isRegular.m\n')
adj = edgeL2adj(canonical_nets(20,'circle'));
assert(isRegular(adj),true)

adj = edgeL2adj(canonical_nets(20,'tree',3));
assert(isRegular(adj),false)

assert(isRegular([0 1; 1 0]),true)
assert(isRegular([0 0; 1 0]),false)
% ================================================

% testing isComplete.m ============================
printf('testing isComplete.m\n')
assert(isComplete([0 1; 1 0]),true)

assert(isComplete(edgeL2adj(directed_cherry)),false)

assert(isComplete(edgeL2adj(undirected_cherry)),false)

randint = randi(10)+10;
adj = ones(randint)-eye(randint);
assert(isComplete(adj),true)
% ================================================

% testing isEulerian.m ===========================
printf('testing isEulerian.m\n')

adj = edgeL2adj(canonical_nets(10,'circle'));
assert(isEulerian(adj),true)

adj = edgeL2adj(canonical_nets(10,'tree',3));
assert(isEulerian(adj),false)
% ================================================

% ================================================
printf('testing isTree.m\n')
adj = edgeL2adj(canonical_nets(randi(10)+10,'tree',2));
assert(isTree(adj),true)

adj = edgeL2adj(canonical_nets(randi(10)+10,'circle'));
assert(isTree(adj),false)
% ================================================

% testing isGraphic.m ============================
printf('testing isGraphic.m\n')
for i=1:100
  adj = giantComponent(random_graph(randi(20)+1,0.5));
  [deg,~,~] = degrees(adj);
  assert(isGraphic(deg) | adj==0)
end
% ================================================

% testing isBipartite.m ==========================
printf('testing isBipartite.m\n')

assert(isBipartite(adj2adjL(bowtie)),false)
assert(isBipartite(edgeL2adjL(undirected_cherry)),true)

even_circle = canonical_nets(2*randi(10),'circle');
assert(isBipartite(edgeL2adjL(even_circle)),true)

odd_circle = canonical_nets(2*randi(10)+1,'circle');
assert(isBipartite(edgeL2adjL(odd_circle)),false)
% ================================================


% testing adj2adjL.m =============================
printf('testing adj2adjL.m\n')
assert(adj2adjL(bowtie),bowtie_adjL')
% ================================================

% testing adjL2adj.m =============================
printf('testing adjL2adj.m\n')

assert(adjL2adj(bowtie_adjL),bowtie)

L = {}; L{1}=[2,3]; L{2}=[]; L{3}=[];
assert(adjL2adj(L),edgeL2adj(directed_cherry))
% ================================================

% testing adj2edgeL.m ============================
printf('testing adj2edgeL.m\n')

assert(sortrows(adj2edgeL(bowtie)),bowtie_edgeL)
assert(adj2edgeL([0 1 1; 0 0 0; 0 0 0]),directed_cherry)
% ================================================

% testing edgeL2adj.m ============================
printf('testing edgeL2adj.m\n')

assert(edgeL2adj(bowtie_edgeL),bowtie)
assert(edgeL2adj(directed_cherry),[0 1 1; 0 0 0; 0 0 0])
% ================================================

% testing adj2inc.m ==============================
printf('testing adj2inc.m\n')

randint = randi(10)+1;
assert(adj2inc(eye(randint)),eye(randint))

assert(adj2inc([0 1 0; 0 1 0; 1 0 0 ]),[-1 0 1; 1 1 0; 0 0 -1])
assert(adj2inc([0 2; 0 0]),[-1 -1; 1 1])  % double edge
% ================================================

% testing inc2adj.m ==============================
printf('testing inc2adj.m\n')

randint = randi(10)+1;
assert(inc2adj(eye(randint))==eye(randint))

adj = ones(3) - eye(3);
assert(inc2adj(adj),adj)

inc = [-1 1; 1 0; 0 -1];  % two edges (1->2, 3->1)
assert(inc2adj(inc)==[0 1 0; 0 0 0; 1 0 0])
% ================================================

% testing adj2str.m ==============================
printf('testing adj2str.m\n')

assert(adj2str(ones(3)-eye(3)),'.2.3,.1.3,.1.2,')
assert(adj2str(eye(3)),'.1,.2,.3,')
assert(adj2str([0 2; 0 0]),'.2,,')
% ================================================

% testing str2adj.m ==============================
printf('testing str2adj.m\n')

assert(ones(3)-eye(3),str2adj('.2.3,.1.3,.1.2,'))
assert(eye(3),str2adj('.1,.2,.3,'))
assert([0 1 0; 0 0 0; 1 0 0 ],str2adj('.2,,.1,'))
% ================================================

% testing adjL2edgeL.m ===========================
printf('testing adjL2edgeL.m\n')

assert(adjL2edgeL({[2,3],[],[]}),directed_cherry)
assert(sortrows(adjL2edgeL(bowtie_adjL)),bowtie_edgeL)
% ================================================

% testing edgeL2adjL.m ===========================
printf('testing edgeL2adjL.m\n')
assert(edgeL2adjL(directed_cherry),{[2,3],[],[]}')
% ================================================

% testing inc2edgeL.m ============================
printf('testing inc2edgeL.m\n')

assert(inc2edgeL([1 0 0; 0 1 0; 0 0 1]),[1 1 1; 2 2 1; 3 3 1])  % three self-loops
assert(inc2edgeL([-1 -1; 1 0; 0 1]),[1 2 1; 1 3 1])
assert(inc2edgeL([-1;1]),[1 2 1])
% ================================================


% testing adj2simple.m ===========================
printf('testing adj2simple.m\n')

assert(adj2simple(rand(6)),ones(6)-eye(6))
assert(adj2simple([0 2 0; 1 0 0; 1 2 0]),[0 1 1; 1 0 1; 1 1 0])
% ================================================


% testing edgeL2simple.m =========================
printf('testing edgeL2simple.m\n')

assert(length(edgeL2simple([1 1 1; 2 2 1; 3 3 1])),0)
assert(sortrows(edgeL2simple([1 2 1; 1 3 2;4 5 1.4])),[1 2 1; 1 3 1; 2 1 1; 3 1 1; 4 5 1; 5 4 1])
% ================================================


% testing symmetrize.m ===========================
printf('testing symmetrize.m\n')
for i=1:20
  adj = random_directed_graph(randi(10)+3,rand);
  assert(isSymmetric(symmetrize(adj)),true)
end

% testing symmetrizeEdgeL.m ======================
fprintf('testing symmetrizeEdgeL.m\n')

for x=1:50
  adj = random_directed_graph(randi(20)+2,rand); % create a random adjacency
  el = adj2edgeL(adj);
  if isempty(el); continue; end
  elsym = symmetrizeEdgeL(el);
  adjsym = edgeL2adj(elsym);
  assert(isSymmetric(adjsym),true)
end
% ================================================


% testing addEdgeWeights.m =======================
fprintf('testing addEdgeWeights.m\n')

assert([1 2 2; 1 3 1; 3 4 3],addEdgeWeights([1 2 1; 1 2 1; 1 3 1; 3 4 2; 3 4 1]))
assert([1 2 2; 2 3 4],addEdgeWeights([1 2 2; 2 3 4]))
% ================================================

% ================================================
printf('testing degrees.m\n')

assert([2 2 3 3 2 2],degrees(bowtie))
assert([2 1 1],degrees(edgeL2adj(directed_cherry)))
assert([2 1 1],degrees(edgeL2adj(undirected_cherry)))

[deg,indeg,outdeg]=degrees(edgeL2adj(directed_cherry));
assert(deg,[2 1 1])
assert(indeg,[0 1 1])
assert(outdeg,[2 0 0])

assert([4 4 4],degrees([0 2 1; 0 0 1; 1 1 0]))
% ================================================


% testing rewire.m ===============================
printf('testing rewire.m\n')

for x=1:100
  
  el = adj2edgeL(random_graph(randi(10)+10,0.4));
  deg = degrees(edgeL2adj(el));
  eln = rewire(el,randi(5));
  degn = degrees(edgeL2adj(eln));
  
  assert(deg,degn)

end
% ================================================


% testing rewireAssort.m =========================
printf('testing rewireAssort.m\n')

for x=1:100
  adj = [0 0; 0 0];
  
  while not(isConnected(adj)); adj = random_graph(randi(10)+10,0.4); end
  el = adj2edgeL(adj);
  eln = rewireAssort(el,randi(5));
  
  assert(pearson(edgeL2adj(eln))>=pearson(edgeL2adj(el))-10^(-7))
end
% ================================================


% testing rewireDisassort.m ======================
printf('testing rewireDisassort.m\n')
for x=1:100
  
  adj = [0 0; 0 0];
  while not(isConnected(adj)); adj = random_graph(randi(10)+10,0.4); end
  el = adj2edgeL(adj);
  eln = rewireDisassort(el,randi(5));

  assert(pearson(edgeL2adj(eln))<=pearson(edgeL2adj(el))+10^(-7))
  
end
% ================================================


% testing aveNeighborDeg.m =======================
printf('testing aveNeighborDeg.m\n')
assert(aveNeighborDeg(undirected_triangle),[2 2 2])
assert(aveNeighborDeg(bowtie),[2.5 2.5 7/3 7/3 2.5 2.5])
% ================================================

% testing sortNodesBySumNeighborDegrees.m ===
printf('testing sortNodesBySumNeighborDegrees.m\n')

assert(sortNodesBySumNeighborDegrees(bowtie),[4,3,6,5,2,1]')  
assert(sortNodesBySumNeighborDegrees([0 1 1; 1 0 0; 1 0 0]),[1, 3, 2]')
% ================================================

% testing sortNodesByMaxNeighborDegree.m ====
printf('testing sortNodesByMaxNeighborDegree.m\n')

assert(sortNodesByMaxNeighborDegree(bowtie),[4,3,6,5,2,1]')
assert(sortNodesByMaxNeighborDegree(edgeL2adj(undirected_cherry)),[1,3,2]')

% testing closeness.m ============================
printf('testing closeness.m\n')
assert(closeness(bowtie)',[1/(1+1+2+3+3), 1/(1+1+2+3+3), 1/(1+1+1+2+2), 1/(1+1+1+2+2), 1/(1+1+2+3+3), 1/(1+1+2+3+3)])
assert(closeness([0 1 1; 1 0 0; 1 0 0]),[0.5 1/3 1/3]')
% ================================================


% testing nodeBetweennessSlow.m ==================
printf('testing nodeBetweennessSlow.m\n')
assert(nodeBetweennessSlow([0 1; 1 0]),[0 0])
assert(nodeBetweennessSlow([1 1; 0 0]),[0 0])
assert(nodeBetweennessSlow([0 1 1; 1 0 0; 1 0 0]),[1/3 0 0])
assert(nodeBetweennessSlow(bowtie),[0 0 0.4 0.4 0 0])
% ================================================

% testing nodeBetweennessFaster.m ================
printf('testing nodeBetweennessFaster.m\n')
assert(nodeBetweennessFaster([0 1; 1 0]),[0 0])
assert(nodeBetweennessFaster([0 1 1; 1 0 0; 1 0 0]),[1/3 0 0])
assert(nodeBetweennessFaster(bowtie),[0 0 0.4 0.4 0 0])

adj = [0 0; 0 0];
for i=1:50
  
  while not(isConnected(adj)); adj = random_graph(randi(10)+5,rand); end
  assert(nodeBetweennessSlow(adj),nodeBetweennessFaster(adj))
  
end
% ================================================


% testing edgeBetweenness.m ======================
printf('testing edgeBetweenness.m\n')

eb_bowtie = adj2edgeL(bowtie);
eb_bowtie(:,3) = [1/30; 4/30; 1/30; 4/30; 4/30; 4/30; 9/30; 9/30; 4/30; 4/30; 4/30; 1/30; 4/30; 1/30];
assert(edgeBetweenness(bowtie),eb_bowtie)

if not(isequal(edgeBetweenness(bowtie),eb_bowtie)); fprintf('edge_betweenness.m does not work with bowtie\n'); end

assert(edgeBetweenness(undirected_triangle),[2 1 1/6; 3 1 1/6; 1 2 1/6; 3 2 1/6; 1 3 1/6; 2 3 1/6])
if not(isequal(edgeBetweenness(undirected_triangle),[2 1 1/6; 3 1 1/6; 1 2 1/6; 3 2 1/6; 1 3 1/6; 2 3 1/6])); fprintf('edge_betweenness.m does not with an undirected_triangle\n'); end

assert(edgeBetweenness([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]),[2 1 1/12; 3 1 1/6; 1 2 1/12; 3 2 1/6; 1 3 1/6; 2 3 1/6; 4 3 3/12; 3 4 3/12])
if not(isequal(edgeBetweenness([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]),[2 1 1/12; 3 1 1/6; 1 2 1/12; 3 2 1/6; 1 3 1/6; 2 3 1/6; 4 3 3/12; 3 4 3/12])); fprintf('edge_betweenness.m does not work with triangle+edge\n'); end
% ================================================


% testing eigenCentrality.m ======================
printf('testing eigenCentrality.m\n')
[v,~]=eig([0 1 1; 1 0 1; 1 1 0]);
assert(eigenCentrality([0 1 1; 1 0 1; 1 1 0]),v(:,3))
% ================================================

% testing clustCoeff.m ===========================
printf('testing clustCoeff.m\n')
assert(clustCoeff(undirected_triangle),1)
assert(clustCoeff(edgeL2adj(undirected_cherry)),0)
assert(clustCoeff(edgeL2adj(canonical_nets(randi(10)+5,'tree',2))),0)
[C1,C2] = clustCoeff(bowtie);
assert([C1,C2],[1/3,7/9])
% ================================================


% testing weightedClustCoeff.m ===================
printf('testing weightedClustCoeff.m\n')
randint = randi(20);
assert(length(weightedClustCoeff(random_graph(randint+5,rand))),randint+5)
% ================================================


% testing pearson.m ==============================
printf('testing pearson.m\n')
assert( pearson( edgeL2adj( canonical_nets(randi(5)+5,'star') ) ) ,-1 )
% ================================================


% testing richClubMetric.m =======================
printf('testing richClubMetric.m\n')
assert(richClubMetric(random_graph(randi(5)+5,rand),12),0)
assert(richClubMetric(bowtie,2),linkDensity(bowtie))

mat = [0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0];
assert(richClubMetric(mat,2),1)
% ================================================


% testing sMetric.m ==============================
fprintf('testing sMetric.m\n')
assert(sMetric(undirected_triangle),2*12)
assert(sMetric(bowtie),2*41)
assert(sMetric(edgeL2adj(directed_cherry)),4)
% ================================================