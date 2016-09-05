% Function to find cycles of size 3.
% Note 1 (EY): Works for directed and undirected graphs.
% Note 2 (EY): Very slow, runs at O(N^3) time with brute force but works
% Note 3 (GB): Reduced time a bit by iterating through neighbor nodes instead of all nodes.
% Note 4 (GB): If there are self-loops, cycles can contain the same node multiple times. 
%
% INPUTs: adjacency matrix, (square matrix of zeros and ones)
% OUTPUTs: list (cell) of 3-cycles;
% 
% Other routines used: adj2adjL.m
% Erdem Yilmaz, January 20, 2016
% GB: Last updated, September 5 2016


function res = loops3rev2(A)


sA = size(A);
if (sA(1) ~= sA(2))
    disp('loops3rev2(): Error - Input must be an adjacency (square) matrix.');
    return;
end

AL = adj2adjL(A);


c=0;
tloops = {};

for i = 1:sA(1)
    for j = 1:length(AL{i})
        j = AL{i}(j);
        for k = 1:length(AL{j})

           k = AL{j}(k);

           if A(k,i)==1

               	du = sort([i j k]);
               	du = strcat( num2str(du(1)),'-',num2str(du(2)),'-',num2str(du(3)) );

                if sum(ismember(tloops,du))==0
                    if (i~=k&&k~=j)
                        c=c+1;
                        tloops{c} = du;    
                    end

                end
                      
           end
        end
    end
end
res = tloops;
if (c==0)
     disp('No cycles of size three.');
else
    disp(strcat('Number of 3-cycles:',num2str(c)));
end
