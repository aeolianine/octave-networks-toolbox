% The eigenvalues of the Laplacian of the graph.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: laplacian eigenvalues, sorted
%
% Other routines used: laplacianMatrix.m
% GB: last updated, Oct 10 2012

function s=graphSpectrum(adj)

[~,D]=eig(laplacianMatrix(adj));
s=-sort(-diag(D)); % sort in decreasing order