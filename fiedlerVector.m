% The vector corresponding to the second 
% smallest eigenvalue of the Laplacian matrix
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: fiedler vector, nx1
%
% Other routines used: laplacianMatrix.m
% GB: last updated, Oct 10 2012

function fv=fiedlerVector(adj)

[V,D]=eig(laplacianMatrix(adj));
[~,Y]=sort(diag(D));
fv=V(:,Y(2));