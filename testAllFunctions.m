% Test code for "Octave routines for Network Analysis"

clear all
close all

% Set of test graphs, in various formats =========
one_double_edge = [0 2; 2 0]; 
bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];   % 'adj'
disconnected_bowtie =[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];   % 'adj'
bowtie_edgeL = [1,2,1; 1,3,1; 2,3,1; 3,4,1; 4,5,1; 4,6,1; 5,6,1]; % 'edgelist'
bowtie_edgeL = sortrows(symmetrize_edgeL(bowtie_edgeL));          % 'edgelist'
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