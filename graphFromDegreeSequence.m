% Constructing a graph from a given degree sequence: deterministic
% Note: This is the Havel-Hakimi algorithm.
%
% Inputs: a graphic degree sequence, [d1,d2, ... dn],
%             where di is the degree of the ith node
% Outputs: adjacency matrix, nxn
%
% Last updated: Oct 21 2012

function adj = graphFromDegreeSequence(seq)

    adj = zeros(length(seq));

    while sum(seq) > 0 % while there are still stubs to connect

        % order stubs by decreasing number of degrees left
        [sorted, I] = sort(-seq);
        n1 = I(1);

        for x = 1:-sorted(1)
            n2 = I(x + 1);
            adj(n1, n2) = adj(n1, n2) + 1;
            adj(n2, n1) = adj(n2, n1) + 1;
            seq(n1) = seq(n1) - 1;
            seq(n2) = seq(n2) - 1;
        end

    end


%!test
%! for x=1:40
%!   adj = [0 1; 0 0];
%!   while not(isConnected(adj)); adj = randomGraph(randi(50)+50,rand); end
%!   adjr=graphFromDegreeSequence(degrees(adj));
%!   assert(isSimple(adjr),true)
%!   assert(degrees(adj),degrees(adjr))
%! end
