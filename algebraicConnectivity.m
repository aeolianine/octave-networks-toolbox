% The algebraic connectivity of a graph:
% the second smallest eigenvalue of the Laplacian
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: algebraic connectivity
%
% Other routines used: graphSpectrum.m
% Last updated: Oct 10 2012

function a = algebraicConnectivity(adj)

    s = graphSpectrum(adj);
    a = s(length(s) - 1);


%!test
%! adj = randomGraph(randi(50)+10,rand);
%! assert(length(algebraicConnectivity(adj)),1)
%! assert( algebraicConnectivity(adj), graphSpectrum(adj)(length(adj)-1) ) 
