% The minimum vertex eccentricity is the graph radius.
%
% Inputs: adjacency matrix (nxn)
% Outputs: graph radius
%
% Other routines used: vertexEccentricity.m
% Last updated: Oct 10 2012

function Rg = graphRadius(adj)

    Rg = min(vertexEccentricity(adj));


%!test
%!shared T, el, adj
%! T = load_test_graphs();
%!assert(graphRadius(T{4}{2}),2)
%!assert(graphRadius(edgeL2adj(T{11}{2})),1)

%! el = canonicalNets(randi(10)+10,'line');
%! adj = edgeL2adj(el);
%!assert(graphRadius(adj),(size(adj,1)-mod(size(adj,1),2))/2)

%!demo
%! bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
%! graphRadius(bowtie)
%! adj = [0 1; 0 0]; 
%! graphRadius(adj)