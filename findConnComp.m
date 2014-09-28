% Algorithm for finding connected components in a graph
% Note: Valid for undirected graphs only.
%
% INPUTS: adj - adjacency matrix, nxn
% OUTPUTS: a list of the components comp{i}=[j1,j2,...jk]
%
% Other routines used: findConnCompI.m, degrees.m
% GB: last updated, September 22, 2012


function comp_mat = findConnComp(adj)

[deg,~,~]=degrees(adj);            % degrees
comp_mat={};                       % initialize components matrix

for i=1:length(deg)
    if deg(i)>0
        done=0;
        for x=1:length(comp_mat)
            if length(find(comp_mat{x}==i))>0   % i in comp_mat(x).mat
                done=1;
                break
            end
        end
        if not(done)
            comp=findConnCompI(adj,i);
            comp_mat{length(comp_mat)+1}=comp;
        end
        
    elseif deg(i)==0
        comp_mat{length(comp_mat)+1}=[i];
    end
    
end