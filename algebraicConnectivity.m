% The algebraic connectivity of a graph: 
% the second smallest eigenvalue of the Laplacian
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: algebraic connectivity
%
% Other routines used: graphSpectrum.m
% GB: last updated, Oct 10 2012

function a=algebraicConnectivity(adj)

s=graphSpectrum(adj);
a=s(length(s)-1);