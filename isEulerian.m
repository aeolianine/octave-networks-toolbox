% Check if a graph is Eulerian, i.e. it has an Eulerian circuit
% "A connected undirected graph is Eulerian if and only if every graph vertex has an even degree."
% "A connected directed graph is Eulerian if and only if every graph vertex has equal in- and out- degree."
% Note 1: Assume that the graph is connected.
% Note 2: If there is an Eulerian trail, it is reported.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: Boolean variable, 0 or 1
%
% Other routines used: degrees.m, isDirected.m
% GB: last updated, Sep 23, 2012

function S=isEulerian(adj)

S=false;

[degs,indeg,outdeg]=degrees(adj);
odd=find(mod(degs,2)==1);

if not(isDirected(adj)) && isempty(odd) % if undirected and all degrees are even
  S=true;

elseif isDirected(adj) && indeg==outdeg % directed and in-degrees equal out-degrees
  S=true;

elseif numel(odd)==2
  printf('isEulerian.m: There is an Eulerian trail from node %2i to node %2i\n',odd(1),odd(2));

end