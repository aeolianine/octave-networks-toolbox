% Constructing edge lists for simple canonical graphs, ex: trees and lattices.
%
% INPUTS: number of nodes, network type, branch factor (for trees only).
%         Network types can be 'line','cycle','star','btree','tree',
%            'hierarchy','trilattice','sqlattice','hexlattice', 'clique'
% OUTPUTS: edge list (mx3)
%
% Note: Produces undirected graphs, i.e. symmetric edge lists.
% Other functions used: symmetrizeEdgeL.m, adj2edgeL.m
% GB: last updated: Oct 27 2012


function el=canonicalNets(n,type,b)

if strcmp(type,'line')
    el=buildLine(n);
    
elseif strcmp(type,'cycle')
    el=buildCycle(n);
    
elseif strcmp(type,'star')
    el=buildStar(n);
    
elseif strcmp(type,'clique')
    el=buildClique(n);
    
elseif strcmp(type,'btree')
    el=buildBinaryTree(n);
    
elseif strcmp(type,'tree')
    el=buildTree(n,b);
    
elseif strcmp(type,'hierarchy')
    el=buildHierarchy(n,b);
    
elseif strcmp(type,'trilattice')
    el=buildTriangularLattice(n);
    
elseif strcmp(type,'sqlattice')
    el=buildSquareLattice(n);
    
elseif strcmp(type,'hexlattice')
    el=buildHexagonalLattice(n);
    
else
  printf('invalid network type; see canonicalNets.m header\n');
  el = 'invalid network type';
  return
end

%---- canonical nets functions -----------------------------------

function el_line=buildLine(n) % line .............................

el_line = [[1:n-1]' [2:n]' ones(n-1,1)];
el_line = symmetrizeEdgeL(el_line);


function el_cyc=buildCycle(n) % cycle ...........................

el_cyc = [[1:n-1]' [2:n]' ones(n-1,1)];
el_cyc = [el_cyc; 1 n 1];
el_cyc = symmetrizeEdgeL(el_cyc);


function el_star=buildStar(n) % star .............................

el_star = [ones(n-1,1) [2:n]' ones(n-1,1)]; 
el_star = symmetrizeEdgeL(el_star);


function el_clique = buildClique(n) % clique .....................

% clique:= complete graph with "n" nodes
el_clique = adj2edgeL(ones(n)-eye(n));


function el_bt=buildBinaryTree(n) % binary tree ..................

el_bt=[];

for i=2:n
  if mod(i,2)==0
    el_bt=[el_bt; i i/2 1];
  else
    el_bt=[el_bt; i (i-1)/2 1];
  end
end

el_bt=symmetrizeEdgeL(el_bt);


function el=buildTree(n,b)  % a general tree .....................

% tree with n "nodes" and branch factor "b"

nodes=1;
queue=1;
el=[];

while nodes < n
    if nodes+b > n
        % add (n-nodes) to first stack member
        for i=1:n-nodes
            el=[el; queue(1) i+nodes 1];
        end
        nodes=n;
    else
        % add b new edges: 
        for bb=1:b
            el=[el; queue(1) nodes+bb 1];
            queue=[queue nodes+bb];
        end
        queue=queue(2:length(queue)); % remove first member
        nodes=nodes+b;                % update node number
    end
end

el=symmetrizeEdgeL(el);


function el=buildHierarchy(n,b)  % a hierarchy ...................

% build a tree with n nodes and b as a branch factor, 
%             where nodes on one level are connected
% INPUTs: number of nodes - n, branch factor - b
% OUTPUTs: edge list


L=ceil(log(1+n*(b-1))/log(b));   % L=log(1+n(b-1))/log(b)
el=buildTree(n,b);               %  build the base

% add the cross-level links
for k=1:L
    
    start_node=1+round((b^(k-1)-1)/(b-1));
  
    for bb=1:b^(k-1)-1
        if start_node+bb-1>n || start_node+bb>n
            el=symmetrizeEdgeL(el);
            return
        end
        el=[el; start_node+bb-1, start_node+bb,1];
    end
    
end

el=symmetrizeEdgeL(el);



function el_tr=buildTriangularLattice(n)  % triangular lattice

% build triangular lattice with n nodes
% as close to a "square" shape as possible

el_tr=[];  % initialize edge list

x=factor(n);

if numel(x)==1
  el_tr=[el_tr; n 1 1];
  n=n-1;
  x=factor(n);
end

if mod(numel(x),2)==0  % if there's an even number of factors, split in two

  f1=prod(x(1:numel(x)/2)); 
  f2=prod(x(numel(x)/2+1:numel(x)));

elseif mod(numel(x),2)==1
  
  f1=prod(x(1:(numel(x)+1)/2)); 
  f2=prod(x((numel(x)+1)/2+1:numel(x)));

end

% the lattice will be f1xf2
% inner mesh
for i=2:f1-1
  for j=2:f2-1
    % (i,j)->f2*(i-1)+j
    el_tr=[el_tr; f2*(i-1)+j f2*(i-2)+j 1];
    el_tr=[el_tr; f2*(i-1)+j f2*(i)+j 1];
    el_tr=[el_tr; f2*(i-1)+j f2*(i-1)+j-1 1];
    el_tr=[el_tr; f2*(i-1)+j f2*(i-1)+j+1 1];
    el_tr=[el_tr; f2*(i-1)+j f2*i+j+1 1]; % added for tri lattice
  end
end

% four corners
el_tr=[el_tr; 1 2 1];
el_tr=[el_tr; 1 f2+1 1];
el_tr=[el_tr; 1 f2+2 1]; % added for tri lattice
el_tr=[el_tr; f2 f2-1 1];
el_tr=[el_tr; f2 2*f2 1];
el_tr=[el_tr; f2*(f1-1)+1 f2*(f1-2)+1 1];
el_tr=[el_tr; f2*(f1-1)+1 f2*(f1-1)+2 1];
el_tr=[el_tr; f1*f2 f2*(f1-1) 1];
el_tr=[el_tr; f1*f2 f2*f1-1 1];

% four walls
for j=2:f2-1
  el_tr=[el_tr; j j-1 1];
  el_tr=[el_tr; j j+1 1];
  el_tr=[el_tr; j f2+j 1];
  el_tr=[el_tr; j f2+j+1 1]; % added for tri lattice
  
  el_tr=[el_tr; f2*(f1-1)+j f2*(f1-1)+j-1 1]; 
  el_tr=[el_tr; f2*(f1-1)+j f2*(f1-1)+j+1 1];
  el_tr=[el_tr; f2*(f1-1)+j f2*(f1-2)+j 1];
end
for i=2:f1-1
  el_tr=[el_tr; f2*(i-1)+1 f2*(i-2)+1 1];
  el_tr=[el_tr; f2*(i-1)+1 f2*i+1 1];
  el_tr=[el_tr; f2*(i-1)+1 f2*(i-1)+2 1];
  el_tr=[el_tr; f2*(i-1)+1 f2*i+2 1]; % added for tri lattice
  
  el_tr=[el_tr; f2*i f2*(i-1) 1];
  el_tr=[el_tr; f2*i f2*(i+1) 1];
  el_tr=[el_tr; f2*i f2*i-1 1];
end

el_tr=symmetrizeEdgeL(el_tr);


function [el_sq]=buildSquareLattice(n)   % square latice

el_sq=[];

x=factor(n);

if numel(x)==1
  el_sq=[el_sq; n 1 1];
  n=n-1;
  x=factor(n);
end

if mod(numel(x),2)==0
  
  f1=prod(x(1:numel(x)/2)); 
  f2=prod(x(numel(x)/2+1:numel(x)));

elseif mod(numel(x),2)==1
  
  f1=prod(x(1:(numel(x)+1)/2)); 
  f2=prod(x((numel(x)+1)/2+1:numel(x)));

end


% the lattice will be f1xf2
% inner mesh
for i=2:f1-1
    for j=2:f2-1
        % (i,j)->f2*(i-1)+j
        el_sq=[el_sq; f2*(i-1)+j f2*(i-2)+j 1];
        el_sq=[el_sq; f2*(i-1)+j f2*(i)+j 1];
        el_sq=[el_sq; f2*(i-1)+j f2*(i-1)+j-1 1];
        el_sq=[el_sq; f2*(i-1)+j f2*(i-1)+j+1 1];
    end
end

% four corners
el_sq=[el_sq; 1 2 1];
el_sq=[el_sq; 1 f2+1 1];
el_sq=[el_sq; f2 f2-1 1];
el_sq=[el_sq; f2 2*f2 1];
el_sq=[el_sq; f2*(f1-1)+1 f2*(f1-2)+1 1];
el_sq=[el_sq; f2*(f1-1)+1 f2*(f1-1)+2 1];
el_sq=[el_sq; f1*f2 f2*(f1-1) 1];
el_sq=[el_sq; f1*f2 f2*f1-1 1];

% four walls
for j=2:f2-1
    el_sq=[el_sq; j j-1 1];
    el_sq=[el_sq; j j+1 1];
    el_sq=[el_sq; j f2+j 1];
  
    el_sq=[el_sq; f2*(f1-1)+j f2*(f1-1)+j-1 1]; 
    el_sq=[el_sq; f2*(f1-1)+j f2*(f1-1)+j+1 1];
    el_sq=[el_sq; f2*(f1-1)+j f2*(f1-2)+j 1];
end
for i=2:f1-1
    el_sq=[el_sq; f2*(i-1)+1 f2*(i-2)+1 1];
    el_sq=[el_sq; f2*(i-1)+1 f2*i+1 1];
    el_sq=[el_sq; f2*(i-1)+1 f2*(i-1)+2 1];
  
    el_sq=[el_sq; f2*i f2*(i-1) 1];
    el_sq=[el_sq; f2*i f2*(i+1) 1];
    el_sq=[el_sq; f2*i f2*i-1 1];
end

el_sq=symmetrizeEdgeL(el_sq);


function el_hex=buildHexagonalLattice(n)  % hexagonal lattice .......


% construct subgraph of the triangular lattice, f1xf2
el_hex=[];
x=factor(n);

if numel(x)==1
  el_hex=[el_hex; n 1 1];
  n=n-1;
  x=factor(n);
end

if mod(numel(x),2)==0
  
  f1=prod(x(1:numel(x)/2));
  f2=prod(x(numel(x)/2+1:numel(x)));
  
elseif mod(numel(x),2)==1
  
  f1=prod(x(1:(numel(x)+1)/2)); 
  f2=prod(x((numel(x)+1)/2+1:numel(x)));
  
end

% the lattice will be f1xf2
% inner mesh
fmax=max(f1,f2);
fmin=min(f1,f2);

for ff=1:fmin  % from 1 to fmin
  
  % rows are from ff*fmax+1 to (ff+1)*fmax - 1
  %           (ff+1)*fmax+1 to (ff+2)*fmax - 2
  for gg=1:fmax % in range(fmax):
    
    if gg<fmax  % if it's not the last node in the row
      el_hex=[el_hex; (ff-1)*fmax+gg,(ff-1)*fmax+gg+1,1];
    end

    if ff<fmin && mod(gg,4)==1      %  gg%4==1
      
      % connect (ff-1)*fmax+gg to ff*fmax+gg+1
      % and (ff-1)*fmax+gg+3 to ff*fmax+gg+2
      n1=(ff-1)*fmax+gg;
      n2=ff*fmax+gg+1;
    
      if n1<fmin*fmax && n2<fmin*fmax && gg+1<=fmax
        el_hex=[el_hex; n1 n2 1];
      end
      
      n1=(ff-1)*fmax+gg+3;
      n2=ff*fmax+gg+2;
      
      if n1<fmin*fmax && n2<fmin*fmax && gg+3<=fmax
        el_hex=[el_hex; n1 n2 1];
      end
      
    end
  end
end

el_hex=symmetrizeEdgeL(el_hex);