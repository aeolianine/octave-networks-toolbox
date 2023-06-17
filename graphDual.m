% Finds the dual of a graph; a dual is the inverted nodes-edges graph
% This is also called the line graph, adjoint graph or the edges adjacency
%
% INPUTs: adjacency (neighbor) list representation of the graph (see adj2adjL.m)
% OUTPUTs: adj (neighbor) list of the corresponding dual graph and
%          cell array of [original] edges, i.e. the new nodes
%
% Note: This routine only works for undirected, simple graphs.
% Last updated: Sep 23 2012

function [dL, edge_array] = graphDual(L)

    dL = {}; % initialize

    for i = 1:length(L)

        for j = 1:length(L{i})

            if i <= L{i}(j)
                dL{length(dL) + 1} = [];
            end

        end

    end

    edge_array = {};

    for i = 1:length(L)

        for j = 1:length(L{i})% add i,L{i}j to list of nodes

            if i <= L{i}(j)
                edge_array{length(edge_array) + 1} = strcat(num2str(i), '-', num2str(L{i}(j)));
            end

        end

        for j = 1:length(L{i})% add i - L{i}j to list of edges

            for k = j + 1:length(L{i})
                edge1 = strcat(num2str(min([i, L{i}(j)])), '-', num2str(max([i, L{i}(j)])));
                edge2 = strcat(num2str(min([i, L{i}(k)])), '-', num2str(max([i, L{i}(k)])));

                ind_edge1 = find(ismember(edge_array, edge1) == 1);
                ind_edge2 = find(ismember(edge_array, edge2) == 1);

                dL{ind_edge1} = unique([dL{ind_edge1}, ind_edge2]);
                dL{ind_edge2} = unique([dL{ind_edge2}, ind_edge1]);
            end

        end

    end


%!test
%!shared T, gd, gdT, L, LT
%! T = load_test_graphs();
%! gd=graphDual(adj2adjL(T{4}{2}));
%! gdT={};
%! gdT{1}=[2,3]; gdT{2}=[1,3,4]; gdT{3}=[1,2,4]; gdT{4}=[2,3,5,6]; gdT{5}=[4,6,7]; gdT{6}=[4,5,7]; gdT{7}=[5,6];
%!assert(gd,gdT)
%! gd=graphDual(adj2adjL(T{13}{2}));
%! gdT={};
%! gdT{1}=[2,3]; gdT{2}=[1,3]; gdT{3}=[1,2];
%!assert(gd,gdT)
%! L={}; LT={}; L{1}=[2]; L{2}=[1]; LT{1}=[];
%!assert(LT,graphDual(L))

%!demo
%! % cycle3 in adjacency matrix format is [0 1 1; 1 0 1; 1 1 0] 
%! % (an example in which the graph and its dual are the same)
%! cycle3 = {[2, 3]; [1, 3]; [1, 2]};
%! [L, edges] = graphDual(cycle3)

%! % undirected 3−node binary tree
%! tree = {[2,3]; [1]; [1]};
%! [L, edges] = graphDual(tree)