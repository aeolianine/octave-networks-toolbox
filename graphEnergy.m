% Graph energy defined as: the sum of the absolute values of the
%      real components of the eigenvalues of the adjacency matrix.
% Source: Gutman, The energy of a graph, Ber. Math. Statist.
%                  Sekt. Forsch-ungszentram Graz. 103 (1978) 1-22.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: graph energy
%
% Last updated: Oct 10 2012

function G = graphEnergy(adj)

    [~, e] = eig(adj); % e are the eigenvalues
    G = sum(abs(real(diag(e))));


%!test
%! adj = randomGraph(randi(50)+10,rand);
%! assert(length(graphEnergy(adj)),1)
%! [~,e]=eig(adj);  % e are the eigenvalues
%! G=sum(abs(real(diag(e))));
%! assert(G, graphEnergy(adj))

%!demo
%! adj = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! graphEnergy(adj)