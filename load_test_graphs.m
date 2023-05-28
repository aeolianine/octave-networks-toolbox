% Set of 19 test graphs, in various formats.
% Graphs are stored in a list structure where list contains the following information per graph (in the same order)
%    {name, if any; graph representation; type of representation; set of nodes; set of edges; number of nodes; number of edges}
%
% Last updated: May 28 2023

function T = load_test_graphs()

    T = {}; % initialize the test graphs

    T{1} = {'one_directed_edge', [0 1; 0 0], 'adjacency', [1 2], [1 2 1], 2, 1};
    T{2} = {'one_undirected_edge', [0 1; 1 0], 'adjacency', [1 2], [1 2 1; 2 1 1], 2, 1};
    T{3} = {'one_double_edge', [0 2; 2 0], 'adjacency', [1 2], [1 2 2; 2 1 2], 2, 2};

    bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
    T{4} = {'bowtie', bowtie, 'adjacency', 1:6, [1 2; 2 1; 2 3; 3 2; 1 3; 3 1; 3 4; 4 3; 4 5; 5 4; 4 6; 6 4; 5 6; 6 5], 6, 7};

    disconnected_bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
    T{5} = {'disconnected_bowtie', disconnected_bowtie, 'adjacency', 1:6, [1 2; 2 1; 2 3; 3 2; 1 3; 3 1; 4 5; 5 4; 4 6; 6 4; 5 6; 6 5], 6, 6};

    directed_bowtie_edgeL = [1, 2, 1; 1, 3, 1; 2, 3, 1; 3, 4, 1; 4, 5, 1; 4, 6, 1; 5, 6, 1];
    T{6} = {'directed_bowtie_edgeL', directed_bowtie_edgeL, 'edgelist', 1:6, directed_bowtie_edgeL, 1:6, 7};

    bowtie_edgeL = sortrows(symmetrizeEdgeL(directed_bowtie_edgeL));
    T{7} = {'bowtie_edgeL', bowtie_edgeL, 'edgelist', 1:6, bowtie_edgeL, 1:6, 7};

    bowtie_edgeL_loop = [bowtie_edgeL; 4 4 1];
    T{8} = {'bowtie_edgeL_loop', bowtie_edgeL_loop, 'edgelist', 1:6, bowtie_edgeL_loop, 1:6, 8};

    bowtie_adjL = {[2, 3], [1, 3], [1, 2, 4], [3, 5, 6], [4, 6], [4, 5]};
    T{9} = {'bowtie_adjL', bowtie_adjL, 'adjlist', 1:6, bowtie_edgeL, 1:6, 7};

    undirected_tree3 = [1, 2, 1; 2, 1, 1; 1, 3, 1; 3, 1, 1];
    T{10} = {'undirected_tree_3nodes', undirected_tree3, 'edgelist', 1:3, undirected_tree3, 3, 2};

    directed_tree3 = [1, 2, 1; 1, 3, 1];
    T{11} = {'directed_tree_3nodes', directed_tree3, 'edgelist', 1:3, directed_tree3, 3, 2};

    directed_tree3_adjL = {}; directed_tree3_adjL{1} = [2, 3]; directed_tree3_adjL{2} = []; directed_tree3_adjL{3} = [];
    T{12} = {'directed_tree_3_adjL', directed_tree3_adjL, 'adjlist', 1:3, directed_tree3, 3, 2};

    undirected_3cycle = [0 1 1; 1 0 1; 1 1 0];
    T{13} = {'undirected_3cycle', undirected_3cycle, 'adjacency', 1:3, [1 2; 2 1; 1 3; 3 1; 2 3; 3 2], 3, 3};

    undirected_3cycle_selfloops = [1 1 1; 1 1 1; 1 1 0];
    T{14} = {'undirected_3cycle_selfloops', undirected_3cycle_selfloops, 'adjacency', 1:3, [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2], 3, 5};

    undirected_3cycle_incidence = [1 1 0; 1 0 1; 0 1 1];
    T{15} = {'undirected_3cycle_incidence', undirected_3cycle_incidence, 'incidence', 1:3, [1 2; 2 1; 1 3; 3 1; 2 3; 3 2], 1:3, 3, 3};

    directed_3cycle = [0 1 0; 0 0 1; 1 0 0];
    T{16} = {'directed_3cycle', directed_3cycle, 'adjacency', 1:3, [1 2; 2 3; 3 1], 3, 3};

    directed_3cycle_adjL = {[2], [3], [1]};
    T{17} = {'directed_3cycle_adjL', directed_3cycle_adjL, 'adjlist', 1:3, [1 2; 2 3; 3 1], 3, 3};

    fourCycle = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
    T{18} = {'4-cycle', fourCycle, 'adjacency', 1:4, [1 2; 1 4; 2 1; 2 3; 3 2; 3 4; 4 1; 4 3], 4, 4};

    star = canonicalNets(5, 'star');
    T{19} = {'5-star', star, 'edgelist', 1:5, [1 2; 1 3; 1 4; 1 5; 2 1; 3 1; 4 1; 5 1], 5, 4};



%!test
%!assert(length(load_test_graphs()), 19)
%!assert(load_test_graphs(){16}{1}, 'directed_3cycle')

%!demo
%! T = load_test_graphs();
%! T{7}  % show the 7th test graph