% Checks whether a matrix is symmetric (has to be square).
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: boolean variable, {0,1}
%
% Last update: Sep 23, 2012

function S = isSymmetric(mat)

    S = false; % default

    if mat == transpose(mat)
        S = true;
    end


%!test
%! for i=1:100
%!   assert(isSymmetric(randomGraph(randi(5)+20,rand)),true)
%!   adj = randomDirectedGraph(randi(5)+20,rand);
%!   assert(not(isSymmetric(adj)) | adj==zeros(size(adj)) | adj==ones(size(adj)))
%! end

%!demo
%! isSymmetric([0 1 1; 1 0 0; 1 0 0])
%! isSymmetric([0 1 1; 0 0 0; 0 0 0])