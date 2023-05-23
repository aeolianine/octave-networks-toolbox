% The eigenvalues of the Laplacian of the graph.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: laplacian eigenvalues, sorted
%
% Other routines used: laplacianMatrix.m
% Last updated: Oct 10 2012

function s = graphSpectrum(adj)

    [~, D] = eig(laplacianMatrix(adj));
    s=-sort(-diag(D)); % sort in decreasing order


%!test
%! adj = randomGraph(randi(50)+10,rand);
%! assert(length(graphSpectrum(adj)),length(adj))
%! L = laplacianMatrix(adj);
%! [~,d] = eig(L);
%! assert( sort(graphSpectrum(adj)), sort(diag(d)) )
