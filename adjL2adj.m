% Convert an adjacency list to an adjacency matrix.
% 
% INPUTS: adjacency list: length n, where L{i_1}=[j_1,j_2,...]
% OUTPUTS: adjacency matrix nxn
% 
% Note: Assume that if node i has no neighbours, then L{i}=[];
% GB: last updated, Sep 25 2012


function adj=adjL2adj(adjL)

adj = zeros(length(adjL));

for i=1:length(adjL)
    for j=1:length(adjL{i}); adj(i,adjL{i}(j))=1; end
end