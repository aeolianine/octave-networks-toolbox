% Weighted clustering coefficient (edge-weights).
% Source: Barrat et al, The architecture of complex weighted networks.
%
% INPUTS: weighted adjacency matrix, nxn
% OUTPUTs: vector of node weighted clustering coefficients, nx1
%
% Other routines used: degrees.m, kneighbors.m
% Last updated: Sep 30 2012

function wC = weightedClustCoeff(adj)

    [deg, ~, ~] = degrees(adj);
    n = size(adj, 1); % number of nodes
    wC = zeros(n, 1); % initialize weighted clust coeff

    for i = 1:n% across all nodes
        neigh = kneighbors(adj, i, 1);

        if length(neigh) < 2
            continue;
        end

        s = 0;

        for ii = 1:length(neigh)

            for jj = 1:length(neigh)

                if adj(neigh(ii), neigh(jj)) > 0
                    s = s + (adj(i, neigh(ii)) + adj(i, neigh(jj))) / 2;
                end

            end

        end

        wC(i) = s / (deg(i) * (length(neigh) - 1));
    end


%!test
%! randint = randi(20);
%! assert(length(weightedClustCoeff(randomGraph(randint + 5, rand))), randint + 5)

%!assert(weightedClustCoeff([0 2 1; 2 0 0; 1 0 0]), [0 0 0]')
%!assert(weightedClustCoeff([0 2 1; 2 0 1; 1 1 0]), [1 1 1]')

%!demo
%! adj = [0 2 1; 2 0 0; 1 0 0];
%! weightedClustCoeff(adj)
%! adj=[0 2 1 1; 2 0 3 0; 1 3 0 0; 1 0 0 0];
%! weightedClustCoeff(adj)


% % ALTERNATIVE  ...........................
% wadj=adj;
% adj=adj>0;
%
% [wdeg,~,~]=degrees(wadj);
% [deg,~,~]=degrees(adj);
% n=size(adj,1); % number of nodes
% wC=zeros(n,1);
%
% for i=1:n
%     if deg(i)<2; continue; end
%
%     s=0;
%     for ii=1:n
%         for jj=1:n
%             s=s+adj(i,ii)*adj(i,jj)*adj(ii,jj)*(wadj(i,ii)+wadj(i,jj))/2;
%         end
%     end
%
%     wC(i)=s/(wdeg(i)*(deg(i)-1));
% end
% ..........................................
