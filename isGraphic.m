% Check whether a sequence of number is graphic, i.e.
%            a graph with this degree sequence exists
% Source: Erdős, P. and Gallai, T. "Graphs with Prescribed Degrees
%          of Vertices" [Hungarian]. Mat. Lapok. 11, 264-274, 1960.
%
% INPUTs: a sequence (vector) of numbers
% OUTPUTs: boolean, true or false
%
% Note: Not generalized to directed graph degree sequences.
% Last updated: Sep 24, 2012

function B = isGraphic(seq)

    if not(isempty(find(seq <= 0))) || mod(sum(seq), 2) == 1
        % there are non-positive degrees or their sum is odd
        B = false;
        return;
    end

    n = length(seq);
    seq=-sort(-seq); % sort in decreasing order

    for k = 1:n - 1
        sum_dk = sum(seq(1:k));
        sum_dk1 = sum(min([k * ones(1, n - k); seq(k + 1:n)]));

        if sum_dk > k * (k - 1) + sum_dk1
            B = false;
            return;
        end

    end

    B = true;

%!test
%! for i=1:20
%!     adj = giantComponent(randomGraph(randi(20)+1,0.5));
%!     [deg,~,~] = degrees(adj);
%!     assert(isGraphic(deg) | adj==0)
%! end

%!assert(isGraphic([0 1]), false)
%!assert(isGraphic([2 1]), false)
%!assert(isGraphic([1 1 4]), false)
%!assert(isGraphic([1 4 4 100]), false)

%!demo
%! isGraphic([2 2 2])
%! isGraphic([2 1 2])
%! isGraphic([2 2 3 3 2 2])
%! isGraphic([2 2 3 3 2 3])