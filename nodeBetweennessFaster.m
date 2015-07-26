% Betweenness centrality measure: number of shortest paths running through a vertex.
% Reference: Ulrik Brandes, "A Faster Algorithm for Betweenness
% Centrality", Journal of Mathematical Sociology 25(2):163-177,(2001)
% Source: http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf
%
% INPUTS: adjacency or distances matrix, nxn
% OUTPUTS: betweeness vector for all vertices (1xn)
%
% Other routines used: kneighbors.m
% GB: July 18 2015

function betw = nodeBetweennessFaster(adj)

betw = zeros(1,size(adj,1));

for s=1:size(adj,1)  % over all nodes
    
    S = [];   % S <- empty stack
    P = {}; P{size(adj,1)} = [];  % P{w} <- empty list, w in V;

    sigma = zeros(1,size(adj,1)); sigma(s) = 1;
    d = (-1)*ones(1,size(adj,1)); d(s) = 0;
    Q = [];  % Q - empty queue
    Q = [Q s];  % enqueue s -> Q
    
    while length(Q)>0  % while Q not empty do
        
        v = Q(1); Q = Q(2:length(Q)); % dequeue v <- Q;
        S = [v S];  % push v -> S
        
        knei = kneighbors(adj,v,1);
        for ii=1:length(knei) % for each neighbor w of v
            w = knei(ii);
            % w found for the first time?
            if d(w)<0
                Q = [Q w];
                d(w) = d(v) + 1;
            end
            % shortest path to w via v?
            if d(w) == d(v) + 1
                sigma(w) = sigma(w) + sigma(v);
                P{w} = [P{w} v];
            end
        end
    end
    
    delta = zeros(1,size(adj,1));
    % S returns vertices in order of non-increasing distance from s
    while length(S)>0   % while S not empty
        w = S(1); S = S(2:length(S));  % pop w <- S

        for ii=1:length(P{w})
            v = P{w}(ii);
            delta(v) = delta(v) + (sigma(v)/sigma(w))*(1+delta(w));
        end
        if w~=s; betw(w) = betw(w) + delta(w); end
    end
end


% this last step is just additional normalization, and is arbitrary 
betw=betw/(2*nchoosek(size(adj,1),2));