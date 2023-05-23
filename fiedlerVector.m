% The vector corresponding to the second
% smallest eigenvalue of the Laplacian matrix
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: fiedler vector, nx1
%
% Other routines used: laplacianMatrix.m
% Last updated: Oct 10 2012

function fv = fiedlerVector(adj)

    [V, D] = eig(laplacianMatrix(adj));
    [~, Y] = sort(diag(D));
    fv = V(:, Y(2));


%!test
%! adj = randomGraph(randi(50)+10,rand);
%! assert(length(fiedlerVector(adj)),length(adj))
%! [V,D]=eig(laplacianMatrix(adj));
%! [~,Y]=sort(diag(D));
%! fv=V(:,Y(2));
%! assert(fv, fiedlerVector(adj))
