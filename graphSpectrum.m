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

%!demo
%! graphSpectrum([0 1; 0 0])
%! adj = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! graphSpectrum(adj)
%! % Notice that the last value of the bowtie graph spectrum is essentially 0