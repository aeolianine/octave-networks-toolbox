% Calculate the number of star motifs of given (subgraph) size.
% Note 1: Easily extendible to return the actual stars as k-tuples of nodes.
% Note 2: Star of size 1 is the trivial case of a single node.
%
% INPUTs: adjacency list {} (1xn), k - star motif size
% OUTPUTs: number of stars with k nodes (k-1 spokes)
%
% Last updated: Oct 5, 2012

function num = numStarMotifs(adjL, k)

    num = 0;

    for i = 1:length(adjL)

        if length(adjL{i}) >= (k - 1);
            num = num + nchoosek(length(adjL{i}), k - 1);
        end

    end

    
% ALTERNATIVE
% function num = numStarMotifs(adj,k)
%
% % INPUTs: adjacency matrix, k - star motif size
% % OUTPUTs: number of stars with k nodes (k-1 spokes)
%
% [deg,~,~]=degrees(adj);
%
% num=0;
%
% for i=1:length(deg)
%     if deg(i)>=(k-1); 
%         num=num+nchoosek(deg(i),k-1); 
%     end
% end

%!test
%!shared T
%! T = load_test_graphs();
%!assert(numStarMotifs(T{9}{2},3),4+6)
%!assert(numStarMotifs(T{9}{2},4),2)
%!assert(numStarMotifs(T{9}{2},5),0)
%!assert(numStarMotifs(adj2adjL(T{13}{2}),3),3)
%!assert(numStarMotifs(adj2adjL(T{13}{2}),2),6)
%!assert(numStarMotifs(T{9}{2},1),6)   % trivial case

%!demo
%! s = numStarMotifs({[2, 3, 4]}, 3)
%! bowtie_adjL = {[2, 3], [1, 3], [1, 2, 4], [3, 5, 6], [4, 6], [4, 5]};
%! s = numStarMotifs(bowtie_adjL, 4)
%! s = numStarMotifs(bowtie_adjL, 5)
