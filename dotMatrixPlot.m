% Draws the matrix as a column/row sorted square dot-matrix pattern.
%
% INPUTs: adj (nxn) - adjacency matrix representation of the graph
% OUTPUTs: plot
%
% Note 1: Change colors and marker types in lines 64, 71, 78, 85 and 92, 99
% Note 2: Easy to add/remove different node orderings to/from the plot
%
% Other routines used: degrees.m, sortNodesByMaxNeighborDegree.m,
%                      eigenCentrality.m, newmanEigenvectorMethod.m,
%                      nodeBetweennessFaster.m, newmanGirvan.m, closeness.m
% GB: last updated, April 14 2015

function [] = dotMatrixPlot(adj)

n = size(adj,1);
%markersize=ceil(n/50); % scale for plotting purposes
markersize=3;
deg = degrees(adj);

Yd = sortNodesByMaxNeighborDegree(adj); % degree centrality
[~, Yb] = sort(nodeBetweennessFaster(adj)); % node betweenness centrality
[~,Yec] = sort(eigenCentrality(adj)); % eigen-centrality
[~,Yc] = sort(closeness(adj)); % closeness

% sort by module
modules = newmanEigenvectorMethod(adj);
% sort modules by length
mL=zeros(1,length(modules));
for i=1:length(modules); mL(i)=length(modules{i}); end
[~,Yms]=sort(mL);

% sort nodes by degree inside modules
Ym=[];
for mm=1:length(modules)
    module=modules{Yms(mm)};
    deg_module=deg(module);
    [~,Yds]=sort(deg_module);
    module_sorted=module(Yds);
    for xx=1:length(module_sorted)
        Ym=[Ym module_sorted(xx)];
    end
end

mods = newmanGirvan(adj,4);   % enter the expected number of communities

% sort modules by length
mL=zeros(1,length(mods));
for i=1:length(mods); mL(i)=length(mods{i}); end
[~,Yms]=sort(mL);
Ymb = [];
for mm=1:length(mods)
    module=mods{Yms(mm)};
    deg_module=deg(module);
    [~,Yds]=sort(deg_module);
    module_sorted=module(Yds);
    for xx=1:length(module_sorted)
        Ymb=[Ymb module_sorted(xx)];
    end
end

set(gcf,'Color',[1 1 1])
subplot(3,2,1)
spy(adj(Yd,Yd),'k.',markersize)
title('ordered by degree','FontWeight','bold')
axis([0 n 0 n]);
set(gca,'YDir','normal')
axis square

subplot(3,2,2)
spy(adj(Yb,Yb),'k.',markersize)
title('ordered by betweenness','FontWeight','bold')
axis([0 n 0 n]);
set(gca,'YDir','normal')
axis square

subplot(3,2,3)
spy(adj(Yec,Yec),'k.',markersize)
title('ordered by eigen-centrality','FontWeight','bold')
axis([0 n 0 n]);
set(gca,'YDir','normal')
axis square

subplot(3,2,4)
spy(adj(Yc,Yc),'k.',markersize)
title('ordered by closeness','FontWeight','bold')
axis([0 n 0 n]);
set(gca,'YDir','normal')
axis square

subplot(3,2,5)
spy(adj(Ym,Ym),'k.',markersize)
title('ordered by module - eigenvector method','FontWeight','bold')
axis([0 n 0 n]);
set(gca,'YDir','normal')
axis square

subplot(3,2,6)
spy(adj(Ymb,Ymb),'k.',markersize)
title('ordered by module - Newman-Girvan','FontWeight','bold')
axis([0 n 0 n]);
set(gca,'YDir','normal')
axis square