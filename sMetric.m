% The sum of products of degrees across all edges.
% Source: "Towards a Theory of Scale-Free Graphs: Definition, Properties, and Implications", by Li, Alderson, Doyle, Willinger
% Note: The total degree is used regardless of whether the graph is directed or not.
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: s-metric
%
% Other routines used: degrees.m
% Last updated: Oct 1 2012

function s = sMetric(adj)

    [deg, ~, ~] = degrees(adj);
    edges = find(adj > 0);

    s = 0;

    for e = 1:length(edges)
        [i, j] = ind2sub([length(adj), length(adj)], edges(e));
        s = s + deg(i) * deg(j);
    end

% ALTERNATIVE ....................
% [deg,~,~]=degrees(adj);
% el=adj2edgeL(adj);
% 
% s=0;
% for e=1:size(el,1)
%     if el(e,1)==el(e,2)
%         s=s+deg(el(e,1))*deg(el(e,2))*el(e,3)*2;  % count self-loops twice
%     else
%         s=s+deg(el(e,1))*deg(el(e,2))*el(e,3);  % multiply by the weight for edges with weights
%     end
% end
% ................................