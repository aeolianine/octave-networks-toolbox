% Simple implementation of breadth-first-search of a graph.
% Returns a breadth-first-search tree.
%
% INPUTs: adjacency list (nxn), "adjL"
%         start node index, "s"
%         end node index, "t"
% OUTPUTs: BFS tree, in adjacency list {} format (directed)
%          (This is the tree that the algorithm walks in
%          search of the target. If the target is not found,
%          this tree is effectively a spanning tree of the
%          entire graph.
%
% GB: last updated, Nov 8 2014

function T=BFS(adjL,s,t)

discovered=[s];          % mark visited
q=[s];                   % add to queue
T=cell(1,length(adjL));   % initialize search path/tree

if s==t
    printf('BFS(): The start and the end node are the same.\n')
    return
end

while not(isempty(q))
    j=q(1); q=q(2:length(q));  % pop the front
    neigh=adjL{j};
    for nn=1:length(neigh)
        if isempty(find(discovered==neigh(nn)))
            T{j}=[T{j}, neigh(nn)];
            if neigh(nn)==t
                return
            end
            discovered=[discovered, neigh(nn)];
            q=[q, neigh(nn)];
        end
    end
end

printf('BFS(): Target node not found.\n')

