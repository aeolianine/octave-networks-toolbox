% Calculating the Pearson coefficient for a degree sequence.
% Source: "Assortative Mixing in Networks", M.E.J. Newman, Phys Rev Let 2002
%
% INPUTs: M - (adjacency) matrix, nxn (square)
% OUTPUTs: r - Pearson coefficient
%
% Other routines used: degrees.m, numEdges.m, adj2inc.m
% See also pearsonW.m
% Last updated: October 1, 2012

function r = pearson(M)

    [degs, ~, ~] = degrees(M); % get the total degree sequence
    m = numEdges(M); % number of edges in M
    inc = adj2inc(M); % get incidence matrix for convience

    % j,k - remaining degrees of adjacent nodes for a given edge
    % sumjk - sum of all products jk
    % sumjplusk - sum of all sums j+k
    % sumj2plusk2 - sum of all sums of squares j^2+k^2

    % compute sumjk, sumjplusk, sumj2plusk2
    sumjk = 0; sumjplusk = 0; sumj2plusk2 = 0;

    for i = 1:m
        [v] = find(inc(:, i) == 1);
        j = degs(v(1)) - 1; k = degs(v(2)) - 1; % remaining degrees of 2 end-nodes
        sumjk = sumjk + j * k;
        sumjplusk = sumjplusk + 0.5 * (j + k);
        sumj2plusk2 = sumj2plusk2 + 0.5 * (j^2 + k^2);
    end

    % Pearson coefficient formula
    r = (sumjk - sumjplusk^2 / m) / (sumj2plusk2 - sumjplusk^2 / m);


%!test
%!shared T
%! T = load_test_graphs();
%!assert(pearson(edgeL2adj(T{10}{2})), -1)
%!assert(pearson(edgeL2adj(T{19}{2})), -1)
%!assert(pearson(edgeL2adj(canonicalNets(randi(5)+5,'star'))), -1)

%!test
%! % test via pearsonW.m (Whitney routine)
%! adj = [1 1; 1 1];
%! for i=1:50
%!   if isComplete(adj); continue; end
%!   adj = randomGraph(randi(20)+3,rand);
%!   assert( abs(pearson(adj)-pearsonW(adj))<10^(-6)  )
%! end

%!demo
%! adj = [0 1 1 1 1; 1 0 0 0 0; 1 0 0 0 0; 1 0 0 0 0; 1 0 0 0 0];
%! pearson(adj)
%! pearsonW(adj)
%! bowtie = [0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! pearson(bowtie)
%! pearsonW(bowtie)