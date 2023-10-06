% Construct a random graph given a degree distribution.
% The function first generates the degree sequence by
%       drawing numbers from the specified distribution.
%
% INPUTs: number of nodes, n; distribution type, string
%         "distribution" can be: 'uniform', 'normal',
%         'binomial', 'exponential' and 'custom'
%         if 'custom', W has to be specified: W is a set
%         of probabilities, 1x(n-1), where W(i) is the
%         probability of a node having degree "i".
%         sum(W) = 1
% OUTPUTs: adjacency matrix of the random graph, nxn
%
% Note: Requires the statistics package (">pkg load statistics")
% Other routines used: isGraphic.m, weightedRandomSample.m,
%                          randomGraphFromDegreeSequence.m
% Last updated, May 20 2023

function adj = randomGraphDegreeDist(n, distribution, W)

    adj = [];
    Nseq = 1; % ensure the while loop start
    % always make sure Nseq is a graphic sequence

    switch distribution

        case 'uniform'

            % Here the choice of maximum degree is n/2. In practice, it can
            % be any number between 1 and (n-1). However, at certain (higher)
            % densities random graphs are much harder to construct.

            while not(isGraphic(Nseq))
                Nseq = randi(round((n - 1) / 4), 1, n);
            end

        case 'normal'

            % the (n-1)/10 scaling is so that the degree sequence is safely within [1,n-1] but it is arbitrary
            while not(isGraphic(Nseq))
                Nseq = ceil((n - 1) * randn(1, n) / 10 + (n - 1) / 4);
            end

        case 'binomial'

            p = 0.5; % default parameter for binomial distribution

            while not(isGraphic(Nseq))
                Nseq = ceil(binornd(n - 1, p, 1, n));
            end

        case 'exponential'

            while not(isGraphic(Nseq))
                Nseq = ceil(exprnd((n - 1) / 10, 1, n));
            end

        case 'custom'

            if nargin < 3
                printf('Need to specify a custom distribution.\n');
                return;
            end

            while not(isGraphic(Nseq))
                Nseq = weightedRandomSample(n, [1:n - 1], W);
            end

        otherwise

            printf('Invalid degree distribution type.\n')
            return

    end

    adj = randomGraphFromDegreeSequence(Nseq);



%!test
%!shared N
%! pkg load statistics
%! N = randi(80)+30;
%! adj = randomGraphDegreeDist(N,'uniform');
%! assert(numNodes(adj),N)
%! assert(isSimple(adj),true)
%! [xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
%! plot(xpdf, ypdf, color = 'k')
%! text(5,0.01,strcat('constructing graphs with different degree distributions, N=  ', num2str(N)))
%! hold off; hold on;

%!test
%! adj = randomGraphDegreeDist(N,'normal');
%! assert(numNodes(adj),N)
%! assert(isSimple(adj),true)
%! [xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
%! plot(xpdf, ypdf, color = 'b')
%! hold off; hold on;
    
%!test
%! adj = randomGraphDegreeDist(N,'binomial');
%! assert(numNodes(adj),N)
%! assert(isSimple(adj),true)
%! [xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
%! plot(xpdf, ypdf, color = 'y')
%! hold off; hold on;

%!test
%! adj = randomGraphDegreeDist(N,'exponential');
%! assert(numNodes(adj),N)
%! assert(isSimple(adj),true)
%! [xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(degrees(adj),'noplot');
%! plot(xpdf, ypdf, color = 'r')
%! legend('uniform','normal','binomial','exponential')
%! hold off;

%!test
%! adj = randomGraphDegreeDist(6,'custom',[1/5 1/5 1/5 1/5 1/5]);
%! assert(numNodes(adj),6)
%! assert(isSimple(adj),true)

%!test
%! adj = randomGraphDegreeDist(N,'custom');
%! assert(isempty(adj),true)

%!test
%! adj = randomGraphDegreeDist(N,'anything here');
%! assert(isempty(adj),true)


%!demo
%! adj = randomGraphDegreeDist(20, 'uniform');
%! numNodes(adj)
%! degrees(adj)
%! randomGraphDegreeDist(3, 'custom', [1/2 1/2])
%! randomGraphDegreeDist(3, 'custom', [1/2 1/2])
