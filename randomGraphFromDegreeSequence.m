% Constructing a random graph based on a given degree sequence.
% Idea source: Molloy M. & Reed, B. (1995) Random Structures and Algorithms 6, 161-179
% 
% INPUTs: a graphic sequence of numbers, 1xn
% OUTPUTs: adjacency matrix of resulting graph, nxn
% 
% Note: The simple version of this algorithm gets stuck about half
%       of the time, so in this implementation the last problematic
%       edge is rewired.
%
% Other routines used: adj2edgeL.m, rewireThisEdge.m, edgeL2adj.m
% GB: last updated, Oct 25 2012


function adj= randomGraphFromDegreeSequence(Nseq)


stubs=Nseq;                % assign degrees to stubs
adj = zeros(length(Nseq)); % initialize adjacency matrix


old_sum = 0;
cnt=0;

while sum(stubs)>0   % while no more stubs are left to connect
      
  if cnt>5                       % if rewiring did not work 5 times
      
    el = adj2edgeL(adj);
    ind = find(stubs>0);
    
    if length(ind) == 1; elr = rewireThisEdge([el; ind(1) ind(1) 1],ind(1),ind(1));  end
    if length(ind) == 2;  elr = rewireThisEdge([el; ind(1) ind(2) 1; ind(2) ind(1) 1],ind(1),ind(2)); end
    
    if length(ind)>2 || isempty(elr)           % restart algorithm
      printf('randomGraphFromDegreeSequence(): restarting ...\n')
      stubs = Nseq;
      adj = zeros(length(Nseq));
      old_sum = 0;
      cnt=0;
    
    else
      adj = edgeL2adj(elr);  % return matrix with last edge rewired
      return
      
    end
 
  end
 
  
  new_sum = sum(stubs);
  
  if old_sum==new_sum; cnt = cnt+1; end       % no new nodes have been connected, counter+1
  if old_sum~=new_sum; cnt=0; end             % new connections, restart count
      
  
  [~,n1] = max(stubs);                % pick the node with highest number of remaining stubs
  
  old_sum = sum(stubs);
    
  ind = find(stubs>0);
  n2 = ind(randi(length(ind)));
    
  if n1==n2; continue; end            % no self-loops
    
  if adj(n1,n2)>0; continue; end      % no double edges
  adj(n1,n2)=1; adj(n2,n1)=1;
  stubs(n1) = stubs(n1) - 1;
  stubs(n2) = stubs(n2) - 1;
    
end