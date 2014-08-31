% Convert an adjacency matrix to a one-line string representation of a graph.
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: string
%
% Note 1: The nomenclature used to construct the string is arbitrary. 
%                            Here we use   .i1.j1.k1,.i2.j2.k2,....
%                            In '.i1.j1.k1,.i2.j2.k2,....', 
%                            "dot" signifies new neighbor, "comma" next node
% Note 2: Edge weights are not reflected in the string representation.
% Example: [0 1 1; 0 0 0; 0 0 0] <=> .2.3,,,
%
% Other routines used: kneighbors.m
% GB: last updated, Sep 25 2012


function str=adj2str(adj)

str='';
n=length(adj);

for i=1:n
    neigh=kneighbors(adj,i,1);
    for k=1:length(neigh); str=strcat(str,'.',num2str(neigh(k))); end
    str=strcat(str,','); % close this node's neighbors list
end