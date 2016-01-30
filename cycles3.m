% Calculate number of loops/cycles of length 3
%
% INPUTs: adj - adjacency matrix, nxn
% OUTPUTs: L3 - number of triangles (loops of length 3)
%
% Note: Valid for an undirected network.
% GB: last updated, Oct 5, 2012

function l3 = cycles3(adj)

l3 = trace(adj^3)/6;   % trace(adj^3)/3!