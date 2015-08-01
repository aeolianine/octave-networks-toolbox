% Return the total number of edges given the adjacency matrix
% INPUTs: adjacency matrix, nxn
% OUTPUTs: m - total number of edges/links
%
% Note: Valid for both directed and undirected, simple or general graph
% Other routines used: selfLoops.m, isSymmetric.m
% GB, last updated Sep 19, 2012

function m = numEdges(adj)

sl=selfLoops(adj); % counting the number of self-loops

if isSymmetric(adj) && sl==0    % undirected simple graph
    m=sum(sum(adj))/2;
    
elseif isSymmetric(adj) && sl>0
    m=(sum(sum(adj))-sl)/2+sl; % counting the self-loops only once

elseif not(isSymmetric(adj))   % directed graph (not necessarily simple)
    m=sum(sum(adj));
    
end