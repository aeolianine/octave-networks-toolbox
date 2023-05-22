% Weighted random sampling.
%
% INPUTs: number of draws from a discrete distribution (n)
%         possible values to pick from, (P)
%         set of normalized weights/probabilities, (W)
% OUTPUTs: s - set of n numbers drawn from P
%              according to the weights in W
%
% Last updated: Oct 31 2012

function s = weightedRandomSample(n, P, W)

    s = [];

    if abs(sum(W) - 1) > 10^(-8);
        printf('The probabilities do not sum up to 1.\n');
        s = 'probabilities do not sum up to 1';
        return;
    end

    % divide the unit interval into |P| segments each with length W_i
    unit = [0, W(1)];

    for w = 2:length(W)
        unit = [unit W(w) + unit(length(unit))];
    end

    % draw a random number in the unit interval uniformly - where does it fall?
    while length(s) < n

        lb = max(find(unit < rand)); % rand is in [unit(lb), unit(lb+1)]
        s = [s P(lb)]; % pick P(lb)

    end


%!test
%!  for x = 1:20
%!      n = randi(10) + 1; % number of draws
%!      P = [1:randi(7) + 10]; % population to draw from
%!      W = rand(1, length(P));
%!      W = W / sum(W); % normalize to create weights that sum to 1
%!      s = weightedRandomSample(n, P, W);
%!      assert(length(s), n)
%!      assert(intersect(P, s), unique(s))% test that all draws exist in the population
%!  end
