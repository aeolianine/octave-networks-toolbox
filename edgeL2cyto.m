% Convert an edge list structure m x [node 1, node 2, link] to 
%    Cytoscape input format (.txt or any text extension works)
%
% Note: In Cytoscape the column separator option is semi-colon ";". 
%             If desired, this is easy to change below in line 18.
%
% INPUTs: edge list - mx3 matrix, m = number of edges, file name text string
% OUTPUTs: text file in Cytoscape format with a semicolon column separator

function [ ]=edgeL2cyto(el,filename)

nodes=unique([el(:,1)', el(:,2)']);
m = size(el,1);     % number of edges

fid = fopen(filename,'wt','native');

for i=1:m
  fprintf(fid,'    %4i ;  %4i ;  %3.6f\n',el(i,1),el(i,2),el(i,3));
end

fclose(fid);
