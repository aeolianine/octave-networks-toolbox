% Create a k-regular graph.
% Note: No solution if k and n are both odd.
%
% INPUTs: n - # nodes, k - degree of each vertex
% OUTPUTs: el - edge list of the k-regular undirected graph
%
% Other routines used: symmetrizeEdgeL.m
% Last updated: Oct 28 2012

function el = kregular(n, k)

    el = [];

    if k > n - 1
        fprintf('a simple graph with n nodes and k>n-1 does not exist\n'); 
        return; 
    end

    if mod(k, 2) == 1 && mod(n, 2) == 1;
        fprintf('no solution for *n* and *k* both odd\n');
        return;
    end

    half_degree = floor(k / 2); % k/2 if k even, else (k-1)/2

    for node = 1:n

        for kk = 1:half_degree

            node_f = mod(node + kk, n);

            if node_f == 0
                node_f = n;
            end

            if not(ismember([node, node_f, 1], el, 'rows'))
                el = [el; node node_f 1];
            end

            node_b = mod(node - kk, n);

            if node_b == 0
                node_b = n;
            end

            if not(ismember([node, node_b, 1], el, 'rows'))
                el = [el; node node_b 1];
            end

        end

    end

    if mod(k, 2) == 1 && mod(n, 2) == 0
        % connect mirror nodes
        for node = 1:n / 2

            node_m = mod(node + n / 2, n);

            if node_m == 0
                node_m = n;
            end

            if not(ismember([node, node_m, 1], el, 'rows'))
                el = [el; node node_m 1];
            end

        end

    end

    el = symmetrizeEdgeL(el);


%!test
%! for x=1:30
%!  
%!   n = randi(20)+5;   % random integer between 6 and 25
%!   k = randi(n-2)+1;  % randon integer between 2 and n-1
%!   if mod(k,2)==1 && mod(n,2)==1; 
%!     % no solution for this case
%!     continue;  
%!   end  
%!   el = kregular(n,k);
%!   adj = edgeL2adj(el);
%!   assert(degrees(adj),k*ones(1,length(adj)))
%! 
%! end


%!demo
%! n = randi(9)+10;    % pick a random n
%! k = randi([2, n-1]); % pick a random k
%! 
%! % no solution for both k and n odd
%! if mod(k,2)==1 && mod(n,2)==1
%!     n=n-1; 
%! end
%! 
%! el = kregular(n, k);
%! adj = edgeL2adj(el);
%! % check that all degrees equal k
%! assert(degrees(adj), k*ones(1,n))