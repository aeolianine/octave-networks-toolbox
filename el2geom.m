% Plot a graph for which nodes have given coordinates. Apply color
%           scheme and thicker lines if edges have varying weights.
%
% INPUTS: extended edge list el[i,:]=[n1 n2 m x1 y1 x2 y2]
% OUTPUTS: geometry plot, higher-weight links are thicker and lighter in color
%
% Note 1: m - edge weight; (x1,y1) are the Euclidean coordinates of n1, (x2,y2) - n2
% Note 2: Easy to change colors and corresponding edge weight coloring
%
% GB: last updated: December 8, 2012

function []=el2geom(el)

set(gcf,'Color',[1 1 1])
map=colormap('winter');   % ................... CHANGE COLORMAP HERE

for i=1:size(el,1)
  
    % plot line between two nodes
    x1=el(i,4); y1=el(i,5);
    x2=el(i,6); y2=el(i,7);

    % edge weights .............  CHANGE COLORMAP-LINE THICKNESS SCALE HERE
    if el(i,3)<2
        color=map(8,:);
        linew=1;
    elseif el(i,3)>=2 && el(i,3)<3
        color=map(2*8,:);
        linew=2;
    elseif el(i,3)>=3 && el(i,3)<4
        color=map(3*8,:);
        linew=3;
    elseif el(i,3)>=4 && el(i,3)<5
        color=map(4*8,:);
        linew=4;
    elseif el(i,3)>=5 && el(i,3)<6
        color=map(5*8,:);
        linew=5;
    elseif el(i,3)>=6 && el(i,3)<7
        color=map(6*8,:);
        linew=6;
    elseif el(i,3)>=7 && el(i,3)<8
        color=map(7*8,:);
        linew=7;
    elseif el(i,3)>=8
        color=map(8*8,:);
        linew=8;
    end
    line([x1 x2],[y1 y2],'LineWidth',linew,'Color',color)

    hold off; hold on;

end
axis equal
axis off