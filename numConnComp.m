% Calculate the number of connected components using the eigenvalues
%                   of the Laplacian - counting the number of zeros
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: positive integer - number of connected components
%
% Other routines used: graphSpectrum.m
% GB: last updated: September 22, 2012

function nc=numConnComp(adj)

s=graphSpectrum(adj);
nc=numel(find(s<10^(-5)));   % zero eigenvalues are sometimes close to zeros numerically