% Support function for tarjan.m
% "Performs a single depth-first search of the graph, finding all
% successors from the node vi, and reporting all strongly connected
% components of that subgraph."
% See: http://en.wikipedia.org/wiki/Tarjan's_strongly_connected_components_algorithm
% 
% INPUTs: start node, vi;
%         graph structure (list), L
%         tarjan.m variables to update: S, ind, v, GSCC
% OUTPUTs: updated tarjan.m variables: S, ind, v, GSCC
% 
% Note: Contains recursion.
% Other routines used: strongConnComp.m
% GB: last updated, Sep 22 2012

function [GSCC,S,ind,v]=strongConnComp(vi,S,ind,v,L,GSCC)


v(vi).index = ind;                % Set the depth index for vi
v(vi).lowlink = ind;
ind = ind + 1;

S = [vi S];                         % Push vi on the stack

for ll=1:length(L{vi})
  vj = L{vi}(ll);                   % Consider successors of vi 
    
  if isempty(v(vj).index)           % Was successor vj visited? 
    
    [GSCC,S,ind,v]=strongConnComp(vj,S,ind,v,L,GSCC);   % Recursion
    v(vi).lowlink = min([v(vi).lowlink, v(vj).lowlink]);
    
  elseif not(isempty(find(S==vj)))            % Is vj on the stack?
    v(vi).lowlink = min([v(vi).lowlink, v(vj).index]);
    
  end
end


if v(vi).lowlink == v(vi).index     % Is v the root of an SCC?
  
  SCC = [vi];
  while 1
    vj = S(1); S = S(2:length(S));
    SCC = [SCC vj];
     
    if vj==vi; SCC = unique(SCC); break; end
  end
  
  GSCC{length(GSCC)+1} = SCC;
  
end