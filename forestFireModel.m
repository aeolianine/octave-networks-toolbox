% Implementation of the forest fire model by Leskovec et al
%                 "Graphs over Time: Densification Laws, Shrinking 
%                             Diameters and Possible Explanations"
%
% Inputs: forward burning probability p in [0,1], 
%         backward burning ratio r, in [0,inf),
%         T - number of nodes
% Outputs: adjacency list of the constructed (directed) graph 
%
% Other routines used: weightedRandomSample.m
% GB: last updated, November 28, 2012


function L = forestFireModel(T,p,r)

if T==1
  L{1} = []; % a single node, no edges
  return
elseif T>=2
  L{1} = [2]; L{2} = [1]; % start with a single edge 
end


for t=3:T
  
  L{t} = [];                % new node arrives
  w = randi(t-1);  % pick a node from 1->t-1 uniformly at random
  queue=[w]; visited = []; 
  
  while not(isempty(queue))

    L{t} = [L{t} w];   % connect t to w
    w = queue(1); visited = [visited w];

    outlinks = L{w}; inlinks = [];
    for ll=1:length(L)
      if sum(find(L{ll}==w))>0; inlinks = [inlinks ll]; end
    end
  
    N = length(unique([outlinks inlinks]));   
    
    x = geornd(1-p,1,1);  % mean p/(1-p) => (1-p) probability
    y = geornd(1-r*p,1,1);  % mean rp/(1-rp) => (1-rp) probability
     
    s_in = weightedRandomSample(y,inlinks,ones(size(inlinks))/length(inlinks));
    s_out = weightedRandomSample(x,outlinks,ones(size(outlinks))/length(outlinks));

    ws = unique([s_in s_out]);
    
    L{t} = unique([L{t} ws]);  % add as outlinks to new node 
    
    for ii=1:length(ws)
      if sum(find(visited==ws(ii)))==0; queue = [queue ws(ii)]; end
    end
    queue = queue(2:length(queue));  % remove w from queue
  end
  
end


for ll=1:length(L); L{ll} = setdiff(L{ll},ll); end  % remove self-loops