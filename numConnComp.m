% Calculate the number of connected components using the eigenvalues
%                   of the Laplacian - counting the number of zeros
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: positive integer - number of connected components
%
% Other routines used: graphSpectrum.m
% Last updated: September 22, 2012

function nc = numConnComp(adj)

    s = graphSpectrum(adj);
    nc = numel(find(s < 10^(-5))); % zero eigenvalues are sometimes close to zeros numerically



%!test
%!shared T
%! T = load_test_graphs();
%!assert(numConnComp(T{5}{2}),2)

%!test
%! randint = randi(51);
%! Adj=zeros(randint*30);
%! for x=1:randint
%!   adj=randomGraph(30,0.5);
%!   Adj(30*(x-1)+1:30*x,30*(x-1)+1:30*x)=adj;
%! end
%! assert(numConnComp(Adj),randint)

%!demo
%! numConnComp([0 0 0; 0 0 0; 0 0 0])
%! numConnComp([0 1 1; 1 0 1; 1 1 0])
%! adj = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! numConnComp(adj)