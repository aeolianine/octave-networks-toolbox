% Convert adjacency matrix to an incidence matrix
% Note: Valid for directed/undirected, simple/not simple graphs
%
% INPUTs: adjacency matrix, nxn
% OUTPUTs: incidence matrix: n x m (number of edges)
%
% Other routines used: isDirected.m
% GB: last updated, Sep 25 2012

function inc = adj2inc(adj)

n=length(adj); % number of nodes

inc=[]; % initialize incidence matrix

if isDirected(adj)

  for i=1:n
    for j=1:n
      
      % handle self-loops
      if i==j; for x=1:adj(i,j); inc=[inc; zeros(1,length(1:i-1)), 1,zeros(1,length(i+1:n))]; end; continue; end
      
      for x=1:adj(i,j) % add multiple edges if any
        if i<j
          inc=[inc; zeros(1,length(1:i-1)),-1,zeros(1,length(i+1:j-1)),1,zeros(1,length(j+1:n))];
        else
          inc=[inc; zeros(1,length(1:j-1)),1,zeros(1,length(j+1:i-1)),-1,zeros(1,length(i+1:n))];
        end
      end
    
    end
  end

else % undirected

  for i=1:n
    for j=i:n
      
      % handle self-loops
      if i==j; for x=1:adj(i,j); inc=[inc; zeros(1,length(1:i-1)), 1,zeros(1,length(i+1:n))]; end; continue; end
      % add multiple edges if any
      for x=1:adj(i,j); inc=[inc; zeros(1,length(1:i-1)),1,zeros(1,length(i+1:j-1)),1,zeros(1,length(j+1:n))]; end
    
    end
  end
  
end

inc=inc';