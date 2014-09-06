% Convert an edge list of a general graph to the edge list of a 
% simple graph (no loops, no double edges, no edge weights, symmetric)
%
% INPUTS: edge list (mx3), m - number of edges
% OUTPUTs: edge list of the corresponding simple graph
%
% Note: Assumes all node pairs [n1,n2,x] occur once; if else see addEdgeWeights.m
% Other routines used: symmetrizeEdgeL.m
% GB: last updated, Sep 25, 2012

function el=edgeL2simple(el)

el(:,3)=ones(size(el,1),1);  % make all edge weights 1

ind=find(not(el(:,1)-el(:,2)==0));  % indices of the "non-self-loops"
el=el(ind,:);

el = symmetrizeEdgeL(el);