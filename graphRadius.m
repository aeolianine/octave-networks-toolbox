% The minimum vertex eccentricity is the graph radius.
%
% Inputs: adjacency matrix (nxn)
% Outputs: graph radius
%
% Other routines used: vertexEccentricity.m
% GB: last updated, Oct 10 2012

function Rg=graphRadius(adj)

Rg=min( vertexEccentricity(adj) );