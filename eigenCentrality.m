% The ith component of the eigenvector corresponding to the greatest 
% eigenvalue gives the centrality score of the ith node in the network.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: eigen(-centrality) vector, nx1
%
% GB: last updated, Sep 29, 2012

function x=eigenCentrality(adj)

[V,D]=eig(adj);
[~,ind]=max(diag(D));
x=V(:,ind);