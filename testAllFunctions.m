% Test code for "Octave tools for Network Analysis"

clear all
close all
  

% Set of test graphs, in various formats

T = {};  % the test graphs structure
         % name, if any; graph; type of representation; set of nodes; set of edges; number of nodes; number of edges

T{1} = {'one_directed_edge', [0 1; 0 0], 'adjacency', [1 2], [1 2 1], 2, 1};
T{2} = {'one_undirected_edge', [0 1; 1 0], 'adjacency', [1 2], [1 2 1; 2 1 1], 2, 1};
T{3} = {'one_double_edge', [0 2; 2 0], 'adjacency', [1 2], [1 2 2; 2 1 2], 2, 2};

bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
T{4} = {'bowtie', bowtie, 'adjacency', 1:6, [1 2; 2 1; 2 3; 3 2; 1 3; 3 1; 3 4; 4 3; 4 5; 5 4; 4 6; 6 4; 5 6; 6 5], 6, 7};

disconnected_bowtie =[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
T{5} = {'disconnected_bowtie', disconnected_bowtie, 'adjacency', 1:6, [1 2; 2 1; 2 3; 3 2; 1 3; 3 1; 4 5; 5 4; 4 6; 6 4; 5 6; 6 5], 6, 6}; 

directed_bowtie_edgeL = [1,2,1; 1,3,1; 2,3,1; 3,4,1; 4,5,1; 4,6,1; 5,6,1];
T{6} = {'directed_bowtie_edgeL', directed_bowtie_edgeL, 'edgelist', 1:6, directed_bowtie_edgeL, 1:6, 7};

bowtie_edgeL = sortrows(symmetrizeEdgeL(directed_bowtie_edgeL));
T{7} = {'bowtie_edgeL', bowtie_edgeL, 'edgelist', 1:6, bowtie_edgeL, 1:6, 7};

bowtie_edgeL_loop = [bowtie_edgeL; 4 4 1];
T{8} = {'bowtie_edgeL_loop', bowtie_edgeL_loop, 'edgelist', 1:6, bowtie_edgeL_loop, 1:6, 8};

bowtie_adjL = {[2,3],[1,3],[1,2,4],[3,5,6],[4,6],[4,5]};
T{9} = {'bowtie_adjL', bowtie_adjL, 'adjlist', 1:6, bowtie_edgeL, 1:6, 7};

undirected_tree3 = [1,2,1; 2,1,1; 1,3,1; 3,1,1];
T{10} = {'undirected_tree_3nodes', undirected_tree3, 'edgelist', 1:3, undirected_tree3, 3, 2};

directed_tree3 = [1,2,1; 1,3,1];
T{11} = {'directed_tree_3nodes', directed_tree3, 'edgelist', 1:3, directed_tree3, 3, 2};

directed_tree3_adjL = {}; directed_tree3_adjL{1}=[2,3]; directed_tree3_adjL{2}=[]; directed_tree3_adjL{3}=[];
T{12} = {'directed_tree_3_adjL', directed_tree3_adjL, 'adjlist', 1:3, directed_tree3, 3, 2};


undirected_3cycle=[0 1 1; 1 0 1; 1 1 0];
T{13} = {'undirected_3cycle', undirected_3cycle, 'adjacency', 1:3, [1 2; 2 1; 1 3; 3 1; 2 3; 3 2], 3, 3};

undirected_3cycle_selfloops = [1 1 1; 1 1 1; 1 1 0];
T{14} = {'undirected_3cycle_selfloops', undirected_3cycle_selfloops, 'adjacency', 1:3, [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2], 3, 5};

undirected_3cycle_incidence = [1 1 0; 1 0 1; 0 1 1];
T{15} = {'undirected_3cycle_incidence', undirected_3cycle_incidence, 'incidence', 1:3, [1 2; 2 1; 1 3; 3 1; 2 3; 3 2], 1:3, 3, 3};

directed_3cycle=[0 1 0; 0 0 1; 1 0 0];
T{16} = {'directed_3cycle', directed_3cycle, 'adjacency', 1:3, [1 2; 2 3; 3 1], 3, 3};

directed_3cycle_adjL={[2], [3], [1]};
T{17} = {'directed_3cycle_adjL', directed_3cycle_adjL, 'adjlist', 1:3, [1 2; 2 3; 3 1], 3, 3};

fourCycle = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
T{18} = {'4-cycle', fourCycle, 'adjacency', 1:4, [1 2; 1 4; 2 1; 2 3; 3 2; 3 4; 4 1; 4 3], 4, 4};

star = canonicalNets(5,'star');
T{19} = {'5-star', star, 'edgelist', 1:5, [1 2; 1 3; 1 4; 1 5; 2 1; 3 1; 4 1; 5 1], 5, 4};



% add another graph test here ....
% ................................................

% ................................................
% ... graph representation functions .............
% ................................................

% Testing adj2adjL.m .............................
printf('testing adj2adjL.m\n')
tic
assert(adj2adjL( T{4}{2} ),T{9}{2}')     % "bowtie" graph
% assert(adj2adjL( edgeL2adj(T{11}{2}) ),T{12}{2})   % directed binary tree
assert(adj2adjL( T{16}{2} ),T{17}{2}')   % directed 3-cycle
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing adjL2adj.m .............................
printf('testing adjL2adj.m\n')
tic
assert(adjL2adj( T{9}{2} ),T{4}{2} )      % "bowtie" graph
assert(adjL2adj( T{17}{2} ),T{16}{2} )    % directed 3-cycle
assert(adjL2adj( T{12}{2} ), edgeL2adj(T{11}{2}) )
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing adj2edgeL.m ............................
printf('testing adj2edgeL.m\n')
tic
for i=1:length(T)
    if not(strcmp( T{i}{3}, 'adjacency' )); continue; end
    edgeL1 = sortrows( adj2edgeL(T{i}{2}) );
    edgeL2 = sortrows( T{i}{5} );
    
    assert(edgeL1(:,1:2), edgeL2(:,1:2))
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing edgeL2adj.m ............................
printf('testing edgeL2adj.m\n')
tic
for i=1:length(T)
    if not(strcmp( T{i}{3}, 'adjacency' )); continue; end
    edgeL = T{i}{5};
    % adding 1s to get the expected edge list dimensions right
    if size(edgeL)(2)==2
        edgeL = [edgeL ones(size(edgeL)(1),1)];
    end
    assert(T{i}{2}, edgeL2adj( edgeL ))
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing adj2inc.m ..............................
printf('testing adj2inc.m\n')
tic
randint = randi(10)+1;
assert(adj2inc(eye(randint)),eye(randint))

assert(adj2inc(T{13}{2}), T{15}{2})   % directed 3-cycle

% 1->2, 2->2, 3->1
assert(adj2inc([0 1 0; 0 1 0; 1 0 0 ]),[-1 0 1; 1 1 0; 0 0 -1])
assert(adj2inc([0 2; 0 0]),[-1 -1; 1 1])  % directed double edge
assert(adj2inc(T{1}{2}), [-1 1]')          % one directed edge
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing inc2adj.m ..............................
printf('testing inc2adj.m\n')
tic
randint = randi(10)+1;
assert(inc2adj(eye(randint))==eye(randint))

adj = ones(3) - eye(3);
assert(inc2adj(adj),adj)

% 1->2, 2->2, 3->1
assert([0 1 0; 0 1 0; 1 0 0 ],inc2adj([-1 0 1; 1 1 0; 0 0 -1]))
assert([0 2; 0 0],inc2adj([-1 -1; 1 1]))  % directed double edge
assert(T{1}{2}, inc2adj([-1 1]'))          % one directed edge

inc = [-1 1; 1 0; 0 -1];  % two edges (1->2, 3->1)
assert(inc2adj(inc)==[0 1 0; 0 0 0; 1 0 0])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing adj2str.m ..............................
printf('testing adj2str.m\n')
tic
assert(adj2str(ones(3)-eye(3)),'.2.3,.1.3,.1.2,')
assert(adj2str(eye(3)),'.1,.2,.3,')
assert(adj2str([0 2; 0 0]),'.2,,')

assert(adj2str(T{4}{2}),'.2.3,.1.3,.1.2.4,.3.5.6,.4.6,.4.5,')
assert(adj2str(T{16}{2}),'.2,.3,.1,')
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing str2adj.m ..............................
printf('testing str2adj.m\n')
tic
assert(ones(3)-eye(3),str2adj('.2.3,.1.3,.1.2,'))
assert(eye(3),str2adj('.1,.2,.3,'))
assert([0 1 0; 0 0 0; 1 0 0 ],str2adj('.2,,.1,'))

assert('.2.3,.1.3,.1.2.4,.3.5.6,.4.6,.4.5,', adj2str(T{4}{2}))
assert('.2,.3,.1,', adj2str(T{16}{2}))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing adjL2edgeL.m ...........................
printf('testing adjL2edgeL.m\n')
tic
assert(adjL2edgeL(T{12}{2}),T{11}{5})           % directed 3-tree
assert(sortrows(adjL2edgeL(T{9}{2}))(1:14,1:2),sortrows(T{4}{5}))   % bowtie graph
assert(sortrows(adjL2edgeL(T{17}{2}))(1:3,1:2),T{16}{5})     % directed 3-cycle
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing edgeL2adjL.m ...........................
printf('testing edgeL2adjL.m\n')
tic
assert(edgeL2adjL(T{11}{5}),T{12}{2}')
assert(edgeL2adjL(sortrows(T{4}{5})),T{9}{2}')
assert(edgeL2adjL(T{16}{5}),T{17}{2}')
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing inc2edgeL.m ............................
printf('testing inc2edgeL.m\n')
tic
assert(inc2edgeL([1 0 0; 0 1 0; 0 0 1]),[1 1 1; 2 2 1; 3 3 1])  % three self-loops
assert(inc2edgeL([-1 -1; 1 0; 0 1]),[1 2 1; 1 3 1])   % tree 3 nodes
assert(inc2edgeL([-1;1]),T{1}{5})                     % one directed edge
assert(inc2edgeL([ 1;1]),T{2}{5})                     % one undirected edge
assert(inc2edgeL([1 1; 1 1]),[1 2 1; 1 2 1; 2 1 1; 2 1 1])     % one double edge
assert(sortrows(inc2edgeL(T{15}{2})(1:length(T{15}{5}),1:2)),sortrows(T{15}{5}))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing adj2simple.m ...........................
printf('testing adj2simple.m\n')
tic
assert(adj2simple(rand(6)),ones(6)-eye(6))
assert(adj2simple([0 2 0; 1 0 0; 1 2 0]),[0 1 1; 1 0 1; 1 1 0])
assert(isSymmetric(adj2simple(rand(7))),true)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing edgeL2simple.m .........................
printf('testing edgeL2simple.m\n')
tic
assert(length(edgeL2simple([1 1 1; 2 2 1; 3 3 1])),0)
assert(sortrows(edgeL2simple([1 2 1; 1 3 2; 4 5 1.4])),[1 2 1; 1 3 1; 2 1 1; 3 1 1; 4 5 1; 5 4 1])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
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

% Testing addEdgeWeights.m .......................
fprintf('testing addEdgeWeights.m\n')
tic
assert([1 2 2; 1 3 1; 3 4 3], addEdgeWeights([1 2 1; 1 2 1; 1 3 1; 3 4 2; 3 4 1]))
assert([1 2 2; 2 3 4], addEdgeWeights([1 2 2; 2 3 4]))
assert([1 2 1; 2 1 1], addEdgeWeights([1 2 1; 2 1 1]))
assert([1 2 1; 2 1 1], addEdgeWeights([1 2 1; 2 1 1]))
assert([1 2 1; 2 1 2], addEdgeWeights([1 2 1; 2 1 1; 2 1 1]))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% ................................................
% ... basic network theory functions .............
% ................................................


% Testing getNodes.m .............................
fprintf('testing getNodes.m\n')
tic
for i=1:length(T); assert( getNodes(T{i}{2},T{i}{3}), T{i}{4} ); end

% randomized test
for i=1:10
    n = randi(100);
    adj = randomDirectedGraph(n);
    assert(getNodes(randomDirectedGraph(n),'adjacency'), 1:n)
    assert(getNodes(randomGraph(n),'adjacency'), 1:n)
end
assert(strcmp(getNodes([],'rgegaerger'),'invalid graph type'))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing getEdges.m ............................
printf('testing getEdges.m\n')
tic
for i=1:length(T)
    edges1 = sortrows( T{i}{5} );
    edges2 = sortrows( getEdges(T{i}{2},T{i}{3}) );
    
    assert( edges1(size(edges1)(1),1:2), edges2(size(edges2)(1),1:2) )
end
assert(strcmp(getEdges([],'rgegfgdfgrger'),'invalid graph type'))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


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
assert(selfLoops(bowtie),0)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing multiEdges.m ..........................
printf('testing multiEdges.m\n')
tic
assert(multiEdges(T{3}{2}),1)
assert(multiEdges([0 2 1; 2 0 1; 1 1 0],1))  % triangle with one double edge
assert(multiEdges([0 0 1; 2 0 0; 0 1 0]),1)  % directed triangle with 1 double edge
assert(multiEdges(randomGraph(randi(15))),0)
assert(multiEdges([0 0 1; 2 0 0; 0 2 0]),2)  % directed triangle with 2 double edges
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing averageDegree.m .......................
printf('testing averageDegree.m\n')
tic
assert(averageDegree(T{2}{2}),1)
assert(averageDegree(T{4}{2}),2+1.0/3)
assert(averageDegree(T{18}{2}),2)
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


% Testing findConnCompI.m ........................
printf('testing findConnCompI.m\n')
tic
assert(findConnCompI(T{5}{2},1),[1,2,3])
assert(findConnCompI(T{5}{2},2),[1,2,3])
assert(findConnCompI(T{5}{2},3),[1,2,3])
assert(findConnCompI(T{5}{2},4),[4,5,6])
assert(findConnCompI(T{5}{2},5),[4,5,6])
assert(findConnCompI(T{5}{2},6),[4,5,6])
assert(findConnCompI([0 1 0; 1 0 0; 0 0 0],1),[1,2])
assert(findConnCompI([0 1 0; 1 0 0; 0 0 0],2),[1,2])
assert(findConnCompI([0 1 0; 1 0 0; 0 0 0],3),[3])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing findConnComp.m ........................
printf('testing findConnComp.m\n')
tic
assert(findConnComp(T{5}{2}),{[1,2,3],[4,5,6]})

clear modules
modules{1}=[0];
randint = randi(21);
Adj = []; adj = [];

% make up a matrix (Adj) of randint disconnected components (adj)
for x=1:randint
  randsecint = randi(5)+5;
  
  % remember the disconnected components in "modules"
  lastnode = modules{length(modules)}(length(modules{length(modules)}));
  modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
  
  % make sure adj is not empty, is connected and the number of nodes is "randsecint"
  while isempty(adj) || not(isConnected(adj)) || not(length(adj)==randsecint); adj=randomGraph(randsecint,0.5); end

  Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
end

modules=modules(2:length(modules));
assert(findConnComp(Adj),modules)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................


% Testing giantComponent.m ......................
printf('testing giantComponent.m\n')
tic
assert(giantComponent([0 1 0; 1 0 0; 0 0 0],3),[0 1; 1 0])

clear modules
modules{1}=[0];
randint = randi(10)+1;
Adj = []; adj = [];

% make up a matrix (Adj) of randint disconnected components (adj)
for x=1:randint
  randsecint = randi(10)+5;
  lastnode = modules{length(modules)}(length(modules{length(modules)}));
  modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
  % make sure adj is not empty, is connected and the number of nodes is "randsecint"
  while isempty(adj) || not(isConnected(adj)) || not(length(adj)==randsecint); adj=randomGraph(randsecint,0.5); end
  Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
end

modules=modules(2:length(modules));
L = [];
for m=1:length(modules); L = [L, length(modules{m})]; end;
[maxL,maxind] = max(L);
[GC, GCnodes] = giantComponent(Adj);
assert(GC, subgraph(Adj,modules{maxind}))
assert(GCnodes, modules{maxind})
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

% Testing graphComplement.m .....................
printf('testing graphComplement.m\n')
tic
mat = [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 0 1 1; 1 1 0 1 0 0; 1 1 1 0 1 0; 1 1 1 0 0 1];
assert(graphComplement(T{4}{2}),mat)
assert(graphComplement(T{13}{2}),eye(3))
assert(graphComplement([0 1 1; 1 0 0; 1 0 0]), [1 0 0; 0 1 1; 0 1 1])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................

% Testing graphDual.m ...........................
printf('testing graphDual.m\n')
tic
gd=graphDual(adj2adjL(T{4}{2}));
gdT={};
gdT{1}=[2,3]; gdT{2}=[1,3,4]; gdT{3}=[1,2,4]; gdT{4}=[2,3,5,6]; gdT{5}=[4,6,7]; gdT{6}=[4,5,7]; gdT{7}=[5,6];
assert(gd,gdT)

gd=graphDual(adj2adjL(T{13}{2}));
gdT={};
gdT{1}=[2,3]; gdT{2}=[1,3]; gdT{3}=[1,2];
assert(gd,gdT)

L={}; LT={}; L{1}=[2]; L{2}=[1]; LT{1}=[];
assert(LT,graphDual(L))
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

% Testing DFS.m .................................
printf('testing DFS.m\n')
tic
allPaths = DFS(T{1}{2}, 1, 2, allPaths = {}, path = [], upperBound = 2);
assert(allPaths, {[1 2]})
allPaths = DFS(T{1}{2}, 2, 1, allPaths = {}, path = [], upperBound = 2);
assert(allPaths, {})

allPaths = DFS(T{2}{2}, 1, 2, allPaths = {}, path = [], upperBound = 2);
assert(allPaths, {[1 2]})
allPaths = DFS(T{2}{2}, 2, 1, allPaths = {}, path = [], upperBound = 2);
assert(allPaths, {[2 1]})

allPaths = DFS(T{4}{2}, 1, 5, allPaths = {}, path = [], upperBound = 1);
assert(allPaths, {})
allPaths = DFS(T{4}{2}, 1, 5, allPaths = {}, path = [], upperBound = 3);
assert(allPaths, {[1 3 4 5]})

allPaths = DFS(T{5}{2}, 1, 5, allPaths = {}, path = [], upperBound = 4);
assert(allPaths, {})
allPaths = DFS(T{5}{2}, 1, 3, allPaths = {}, path = [], upperBound = 3);
assert(allPaths, {[1 2 3], [1 3]})

allPaths = DFS(T{5}{2}, 1, 3, allPaths = {}, path = [], upperBound = 1);
assert(allPaths, {[1 3]})

allPaths = DFS([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0], 1, 3, allPaths = {}, path = [], upperBound = 3);
assert(allPaths, {[1 2 3], [1 4 3]})
allPaths = DFS([0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0], 4, 2, allPaths = {}, path = [], upperBound = 3);
assert(allPaths, {[4 1 2], [4 3 2]})
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ...............................................

% Testing BFS.m .................................
printf('testing BFS.m\n')
tic
adjL = {1:2, 2:[]};
tt = BFS(adjL, 1, 1);
assert(tt, {[], []})
tt = BFS(adjL, 1, 2);
assert(tt{1}, 2)
assert(length(tt),2)
assert(class(tt),'cell')

tt = BFS(adjL, 2, 2);
assert(tt, {[],[]})
tt = BFS(adjL, 2, 3);
assert(tt, {[],[]})
assert(length(tt),2)
assert(class(tt),'cell')

tt = BFS(adjL, 1, 3);
assert(tt{1}, 2)
assert(tt{2}, [])
assert(length(tt),2)
assert(class(tt),'cell')

tt = BFS(T{9}{2}, 1, 4);
assert(tt{1},[2 3])
assert(tt{2},[])
assert(tt{3},[4])
assert(tt{4},[])
assert(tt{5},[])
assert(tt{6},[])

tt = BFS(T{9}{2}, 2, 6);
assert(tt{2},[1 3])
assert(tt{1},[])
assert(tt{3},[4])
assert(tt{4},[5 6])
assert(tt{5},[])
assert(tt{6},[])

tt = BFS(T{9}{2}, 5, 2);
assert(tt{5},[4 6])
assert(tt{6},[])
assert(tt{4},[3])
assert(tt{3},[1 2])
assert(tt{2},[])
assert(tt{1},[])

tt = BFS(T{9}{2}, 5, 10);
assert(tt{5},[4 6])
assert(tt{6},[])
assert(tt{4},[3])
assert(tt{3},[1 2])
assert(tt{2},[])
assert(tt{1},[])
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
assert(isBipartite(adj2adjL(T{4}{2})),false)
[isit, A, B] = isBipartite(edgeL2adjL(T{10}{2}));
assert(isit,true)
assert(A, 1)
assert(B, [2 3])

% test using the signless Laplacian
[~,E] = eig(signlessLaplacian(edgeL2adj(T{10}{2})));
isit1 = false;
if sum( abs(diag(E))<10^(-10) )==1; isit1 = true; end
assert(isit, isit1)

even_cycle = canonicalNets(2*randi(10),'cycle');
[isit, A, B] = isBipartite(edgeL2adjL(even_cycle));
assert(isit,true)
assert(length(A), length(B))
assert(mod(A,2), ones(1,length(A)))
assert(mod(B,2), zeros(1,length(B)))

% test using the signless Laplacian
[~,E] = eig(signlessLaplacian(edgeL2adj(even_cycle)));
isit1 = false;
if sum( abs(diag(E))<10^(-10) )==1; isit1 = true; end
assert(isit, isit1)

odd_cycle = canonicalNets(2*randi(10)+1,'cycle');
assert(isBipartite(edgeL2adjL(odd_cycle)),false)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)

% test using the signless Laplacian
[~,E] = eig(signlessLaplacian(edgeL2adj(odd_cycle)));
isit1 = false;
if sum( abs(diag(E))<10^(-10) )==1; isit1 = true; end
assert(isBipartite(edgeL2adjL(odd_cycle)), isit1)

% ................................................

% ................................................
% ........ centrality measures ...................
% ................................................


% Testing degrees.m ..............................
printf('testing degrees.m\n')
tic
assert(degrees(T{1}{2}), [1 1])
assert(degrees(T{2}{2}), [1 1])
assert(degrees(T{3}{2}), [2 2])
assert([2 2 3 3 2 2],degrees(T{4}{2}))
assert([2 1 1],degrees(edgeL2adj(T{10}{2})))

[deg,indeg,outdeg]=degrees(edgeL2adj(T{11}{2}));
assert(deg,[2 1 1])
assert(indeg,[0 1 1])
assert(outdeg,[2 0 0])

assert(degrees(T{13}{2}, [2 2 2]))
assert(degrees(T{14}{2}, [3 3 2]))
assert(degrees(T{18}{2}, [3 3 3 3]))

assert([4 4 4],degrees([0 2 1; 0 0 1; 1 1 0]))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
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

% Testing aveNeighborDeg.m .......................
printf('testing aveNeighborDeg.m\n')
tic
assert(aveNeighborDeg(T{13}{2}),[2 2 2])
assert(aveNeighborDeg(T{4}{2}),[2.5 2.5 7/3 7/3 2.5 2.5])
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


% Testing closeness.m ............................
printf('testing closeness.m\n')
tic
assert(closeness(T{4}{2})',[1/(1+1+2+3+3), 1/(1+1+2+3+3), 1/(1+1+1+2+2), 1/(1+1+1+2+2), 1/(1+1+2+3+3), 1/(1+1+2+3+3)])
assert(closeness([0 1 1; 1 0 0; 1 0 0]),[0.5 1/3 1/3]')
assert(closeness(T{13}{2}),[1/(1+1), 1/(1+1), 1/(1+1)]')
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
eb_bowtie = adj2edgeL(T{4}{2});
eb_bowtie(:,3) = [1/30; 4/30; 1/30; 4/30; 4/30; 4/30; 9/30; 9/30; 4/30; 4/30; 4/30; 1/30; 4/30; 1/30];

assert(edgeBetweenness(T{4}{2}),eb_bowtie)
assert(edgeBetweenness(T{13}{2}),[2 1 1/6; 3 1 1/6; 1 2 1/6; 3 2 1/6; 1 3 1/6; 2 3 1/6])
assert(edgeBetweenness([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]),[2 1 1/12; 3 1 1/6; 1 2 1/12; 3 2 1/6; 1 3 1/6; 2 3 1/6; 4 3 3/12; 3 4 3/12])
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing eigenCentrality.m ......................
printf('testing eigenCentrality.m\n')
tic
[v,~]=eig([0 1 1; 1 0 1; 1 1 0]);
assert(eigenCentrality([0 1 1; 1 0 1; 1 1 0]),v(:,3))   % "3" is the number of nodes

[v,~]=eig(T{4}{2});
assert( eigenCentrality( T{4}{2} ),v(:,size(T{4}{2},1)) )

[v,~]=eig(T{13}{2});
assert( eigenCentrality( T{13}{2} ),v(:,size(T{13}{2},1)) )

[v,~]=eig(T{18}{2});
ec = v(:,size(T{18}{2},1));
assert( eigenCentrality( T{18}{2} ), ec )
assert( norm( ec(1)*ones(length(ec),1) - ec) < 1*e^(-20) )

adj = edgeL2adj( canonicalNets(randi(10)+2, 'cycle') );
[v,~]=eig(adj);
ec = v(:,size(adj,1));
assert( eigenCentrality( adj ), ec )
assert( norm( ec(1)*ones(length(ec),1) - ec) < 1*e^(-20) )
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing clustCoeff.m ...........................
printf('testing clustCoeff.m\n')
tic
assert(clustCoeff(T{13}{2}),1)
assert(clustCoeff(edgeL2adj(T{10}{2})),0)
assert(clustCoeff(edgeL2adj(canonicalNets(randi(10)+5,'tree',2))),0)
[aveC,C] = clustCoeff(T{4}{2});
assert(aveC,(4+2/3)/6)
assert(C, [1 1 1/3 1/3 1 1]')
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
  assert( abs(pearson(adj)-pearsonW(adj))<10**(-6)  )
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
[d,p]=dijkstra(T{4}{2},1,5);
assert(d,3)
assert(p,[1,3,4,5])

[d,p]=dijkstra(T{13}{2},3,[]);
assert(d,[1,1,0])
assert(p,{[3,1],[3,2],[3]})

[d,p] = dijkstra(T{18}{2},3,[]);
assert(d,[2,1,0,1]);
assert(p,{[3,2,1],[3,2],[3],[3,4]})
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing shortestPathDP.m .......................
printf('testing shortestPathDP.m\n')
tic
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


% Testing findAllShortestPaths.m .................
printf('testing findAllShortestPaths.m\n')
tic
adjL = {}; adjL{1} = [2]; adjL{2} = [];  % 1-2 edge
[allPaths, shortestPathLength] = findAllShortestPaths(adjL,1,2, 2, allPaths={});
assert(shortestPathLength,1)
assert(allPaths{1},'-1-2')
assert(length(allPaths),1)

[allPaths, shortestPathLength] = findAllShortestPaths(adjL,2,1, 2, allPaths={});
assert(shortestPathLength,2)  % not updated, equal to length(adjL)-1
assert(allPaths, {})
assert(length(allPaths),0)

% two alternative paths (1-2-4 and 1-3-4)
adjL = {}; adjL{1} = [2,3]; adjL{2} = [4]; adjL{3} = [4]; adjL{4}=[];
[allPaths, shortestPathLength] = findAllShortestPaths(adjL,1,4, 5, allPaths={});
assert(shortestPathLength,2)
assert(length(allPaths),2)
assert(allPaths{1},'-1-2-4')
assert(allPaths{2},'-1-3-4')

[allPaths, shortestPathLength] = findAllShortestPaths(adjL,4,2, 4, allPaths={});
assert(shortestPathLength,4) % not updated, equal to length(adjL)-1
assert(length(allPaths),0)

% a one-directional cycle
adjL={}; adjL{1}=[2]; adjL{2}=[3]; adjL{3}=[4]; adjL{4}=[1];
[allPaths, shortestPathLength] = findAllShortestPaths(adjL,1,4, 4, allPaths={});
assert(shortestPathLength,3)
assert(length(allPaths),1)
assert(allPaths{1},'-1-2-3-4')

[allPaths, shortestPathLength] = findAllShortestPaths(adjL,4,3, 4, allPaths={});
assert(shortestPathLength,3)
assert(length(allPaths),1)
assert(allPaths{1},'-4-1-2-3')

% undirected cycle
adjL={}; adjL{1}=[2,4]; adjL{2}=[3,1]; adjL{3}=[2,4]; adjL{4}=[1,3];
[allPaths, shortestPathLength] = findAllShortestPaths(adjL,2,4, 4, allPaths={});
assert(shortestPathLength,2)
assert(length(allPaths),2)
assert(allPaths{1},'-2-1-4')
assert(allPaths{2},'-2-3-4')

[allPaths, shortestPathLength] = findAllShortestPaths(adjL,3,1, 3, allPaths={});
assert(shortestPathLength,2)
assert(length(allPaths),2)
assert(allPaths{1},'-3-2-1')
assert(allPaths{2},'-3-4-1')
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

% Testing diameter.m .............................
printf('testing diameter.m\n')
tic
assert(diameter(T{13}{2}),1)
assert(diameter(T{4}{2}),3)

el=canonicalNets(randi(10)+5,'line');
adj = edgeL2adj(el);
assert(diameter(adj),length(adj)-1)

el=canonicalNets(randi(10)+5,'cycle');
adj = edgeL2adj(el);
assert(diameter(adj),floor(length(adj)/2))
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

% Testing graphRadius.m ..........................
printf('testing graphRadius.m\n')
tic
assert(graphRadius(T{4}{2}),2)
assert(graphRadius(edgeL2adj(T{11}{2})),1)

el = canonicalNets(randi(10)+10,'line');
adj = edgeL2adj(el);
assert(graphRadius(adj),(size(adj,1)-mod(size(adj,1),2))/2)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing distanceDistribution.m .................
printf('testing distanceDistribution.m\n')
tic
assert(distanceDistribution(T{4}{2}),[7/15, 4/15, 4/15, 0, 0])
assert(distanceDistribution(T{13}{2}),[1, 0])
assert(distanceDistribution(edgeL2adj(T{10}{2})),[2/3,1/3])
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

% Testing cycles3.m ...............................
printf('testing cycles3.m\n')
tic
assert(cycles3(T{4}{2}),2)
assert(cycles3(T{18}{2}),0)
assert(cycles3(T{13}{2}),1)
assert(cycles3(edgeL2adj(canonicalNets(randi(10)+3,'btree'))),0)
assert(cycles3(edgeL2adj(canonicalNets(4,'trilattice'))),2)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing loops3rev2.m ...............................
printf('testing loops3rev2.m\n')
tic
assert(loops3rev2(T{4}{2}),{'1-2-3','4-5-6'})
assert(loops3rev2(T{18}{2}),{})
assert(loops3rev2(T{13}{2}),{'1-2-3'})
assert(loops3rev2(edgeL2adj(canonicalNets(randi(10)+3,'btree'))),{})
assert(loops3rev2(edgeL2adj(canonicalNets(4,'trilattice'))),{'1-2-4','1-3-4'})
assert(loops3rev2(T{16}{2}),{'1-2-3'})
assert(loops3rev2(fourCycle),{})


printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% Testing cycle4nodes.m ...............................
printf('testing cycle4nodes.m\n')
tic
assert(cycle4nodes(fourCycle), {'1-2-3-4'})
assert(cycle4nodes(T{4}{2}),{})
c4 = ones(4)-eye(4); % clique of size 4
assert(cycle4nodes(c4),{'1-2-3-4'})
c6 = ones(6)-eye(6); % clique of size 6
assert(length(cycle4nodes(c6)),nchoosek(6,4))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing cycles4.m ...........................
printf('testing cycles4.m\n')
tic
assert(cycles4(fourCycle), 1)
assert(cycles4(T{4}{2}),0)
assert(cycles4(c4),3)
c6 = ones(6)-eye(6); % clique of size 6
printf('---Time ellapsed: %3f in minutes.\n', toc/60)

# ................................................

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


% Testing graphSpectrum.m ........................
printf('testing graphSpectrum.m\n')
tic
adj = randomGraph(randi(50)+10,rand);
assert(length(graphSpectrum(adj)),length(adj))
L = laplacianMatrix(adj);
[~,d] = eig(L);
assert( sort(graphSpectrum(adj)), sort(diag(d)) )
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing algebraicConnectivity.m ................
printf('testing algebraicConnectivity.m\n')
tic
adj = randomGraph(randi(50)+10,rand);
assert(length(algebraicConnectivity(adj)),1)
assert( algebraicConnectivity(adj), graphSpectrum(adj)(length(adj)-1) ) 
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing fiedlerVector.m ........................
printf('testing fiedlerVector.m\n')
tic
adj = randomGraph(randi(50)+10,rand);
assert(length(fiedlerVector(adj)),length(adj))
[V,D]=eig(laplacianMatrix(adj));
[~,Y]=sort(diag(D));
fv=V(:,Y(2));
assert(fv, fiedlerVector(adj))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing graphEnergy.m ..........................
printf('testing graphEnergy.m\n')
tic
adj = randomGraph(randi(50)+10,rand);
assert(length(graphEnergy(adj)),1)
[~,e]=eig(adj);  % e are the eigenvalues
G=sum(abs(real(diag(e))));
assert(G, graphEnergy(adj))
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% ................................................
% ............ auxiliary .........................
% ................................................


% Testing weightedRandomSample.m .................
printf('testing weightedRandomSample.m\n');
tic
for x=1:20
  n=randi(10)+1;  % number of draws
  P = [1:randi(7)+10];  % population to draw from
  W = rand(1,length(P));
  W = W/sum(W);         % normalize to create weights that sum to 1
  s = weightedRandomSample(n,P,W);
  assert(length(s),n)
  assert(intersect(P,s),unique(s))   % test that all draws exist in the population                                     
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................


% ................................................
% ......... graph models .........................
% ................................................


% Testing canonicalNets.m ........................
printf('testing canonicalNets.m\n');
tic
for x=1:10
  N = randi(50)+10;
  
  assert(strcmp(canonicalNets([],'rgegfgdfgrger'),'invalid network type'))
  
  elLine = canonicalNets(N,'line'); % test line
  adj = edgeL2adj(elLine);
  assert(numNodes(adj),N);
  assert(isConnected(adj),true)
  assert(degrees(adj),[1 2*ones(1,N-2) 1])
  assert(isTree(adj),true)
  
  elC = canonicalNets(N,'cycle'); % test cycle
  adj = edgeL2adj(elC);
  assert(numNodes(adj),N);
  assert(isConnected(adj),true)
  assert(degrees(adj),2*ones(1,N))
  assert(isTree(adj),false)
  
  elS = canonicalNets(N,'star'); % test star
  adj = edgeL2adj(elS);
  assert(numNodes(adj),N);
  assert(isConnected(adj),true)
  assert(degrees(adj),[N-1 ones(1,N-1)])
  assert(isTree(adj),true)
    
  elCl = canonicalNets(N,'clique'); % test clique
  adj = edgeL2adj(elCl);
  assert(numNodes(adj),N)
  assert(isComplete(adj),true)
  assert(degrees(adj),(N-1)*ones(1,N))
  
  
  elBT = canonicalNets(N,'btree'); % test binary tree
  adj = edgeL2adj(elBT);
  assert(numNodes(adj),N)
  assert(isTree(adj),true)
  deg = degrees(adj);
  assert(max(deg),3)
  assert(deg(1),2)
  assert(deg(N),1)
  
  b = randi([3,6]);
  elT = canonicalNets(N,'tree',b); % test general tree
  adj = edgeL2adj(elT);
  assert(numNodes(adj),N)
  assert(isTree(adj),true)
  deg = degrees(adj);
  assert(max(deg)<=b+1,true)
  assert(deg(1),b)
  assert(deg(N),1)
  
  assert(edgeL2adj(canonicalNets(N,'tree',2)),edgeL2adj(canonicalNets(N,'btree')))
  
  b = randi([3,6]);
  elH = canonicalNets(N,'hierarchy',b); % test hierarchy
  adj = edgeL2adj(elH);  
  assert(numNodes(adj),N)
  assert(isTree(adj),false)
  deg = degrees(adj);
  assert(max(deg)<=b+3,true)
  assert(deg(1),b)
  
  eltr = canonicalNets(N,'trilattice');  % test triangular lattice
  adj = edgeL2adj(eltr);
  assert(numNodes(adj),N)
  assert(isTree(adj),false)
  assert(isComplete(adj),false)
  assert(max(degrees(adj))<=6,true)
  assert(isConnected(adj),true)

  elsq = canonicalNets(N,'sqlattice');  % test square lattice
  adj = edgeL2adj(elsq);
  assert(numNodes(adj),N)
  assert(isTree(adj),false)
  assert(isComplete(adj),false)
  assert(max(degrees(adj))<=4,true)
  assert(isConnected(adj),true)

  elhex = canonicalNets(N,'hexlattice');  % test hexagonal lattice
  adj = edgeL2adj(elhex);
  assert(numNodes(adj),N)
  assert(isTree(adj),false)
  assert(isComplete(adj),false)
  assert(max(degrees(adj))<=3,true)
  assert(isConnected(adj),true)

  
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
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

% Testing graphFromDegreeSequence.m ..............
printf('testing graphFromDegreeSequence.m\n')
tic
for x=1:40
  adj = [0 1; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(50)+50,rand); end
  adjr=graphFromDegreeSequence(degrees(adj));
  assert(isSimple(adjr),true)
  assert(degrees(adj),degrees(adjr))
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

% Testing buildSmaxGraph.m .......................
printf('testing buildSmaxGraph.m\n');
tic
for x=1:50
  adj  = [];
  while not(isConnected(adj)); adj = randomGraph(20,0.1); end
  sm = sMetric(adj);
  elmax1 = buildSmaxGraph(degrees(adj));
  adjmax1 = symmetrize(edgeL2adj(elmax1));
  
  smax = sMetric(adjmax1);
  assert(degrees(adjmax1),degrees(adj))
  assert(smax>=sm,true)
  
  elmax2 = buildSmaxGraph(degrees(adj));
  assert(elmax2,elmax1)
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

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

% Testing exponentialGrowthModel.m ...............
printf('testing exponentialGrowthModel.m\n')
tic
for x=1:10
  el = exponentialGrowthModel(randi(100));
  adj=edgeL2adj(el);
  assert(isConnected(adj),true)
  assert(isTree(adj),true)
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


% Testing fabrikantModel.m .......................
printf('testing fabrikantModel.m\n')
tic
for x=1:20
  adj = fabrikantModel(randi(30)+10,rand*10,'off');
  assert(isConnected(adj),true)
  assert(isTree(adj),true)
  assert(isSimple(adj),true)
end
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................

% Testing DoddsWattsSabel.m ......................
printf('testing DoddsWattsSabel.m\n')
tic
for x=1:40
  randint = randi(50)+2;
  m = randi(round(randint/4));
  adj = DoddsWattsSabel(randint,2,m,10*rand,10*rand);
  assert(numEdges(adj),m+randint-1)
  assert(isTree(adj),false)
  assert(isConnected(adj),true)
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

% Testing forestFireModel.m ......................
printf('testing forestFireModel.m\n');
tic
for x=1:20
  randint = randi(20)+5;
  L = forestFireModel(randint,rand,10*rand);
  adj = symmetrize(adjL2adj(L));
  assert(isSimple(adj),true)
  assert(randint,numNodes(L))

end
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
[gH,Q]=newmanCommFast(bowtie);
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

% Testing el2geom.m ..............................
fprintf('testing el2geom.m\n')
tic
[el,p] = newmanGastner(1000,0.5,[],'off');
elnew = [];
for e=1:size(el,1); elnew = [elnew; el(e,1), el(e,2), randi(9), p(el(e,1),1), p(el(e,1),2), p(el(e,2),1), p(el(e,2),2)]; end
figure
el2geom(elnew)
printf('---Time ellapsed: %3f in minutes.\n', toc/60)
% ................................................
