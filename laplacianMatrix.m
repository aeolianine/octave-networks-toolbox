% The Laplacian matrix defined for a *simple* graph 
% Def: the difference b/w the diagonal degree and the adjacency matrices
% Note: This is not the normalized Laplacian
%
% INPUTS: adjacency matrix, nxn
% OUTPUTs: Laplacian matrix, nxn
% 
% GB: last updated, Oct 10 2012

function L=laplacianMatrix(adj)

L=diag(sum(adj))-adj;


% NORMALIZED Laplacian .............

% def normLaplacianMatrix(adj):

% n=length(adj);
% deg = sum(adj); % for other than simple graphs, 
%                     use [deg,~,~]=degrees(adj);

% L=zeros(n);
% edges=find(adj>0);
% 
% for e=1:length(edges)
%     [ii,jj]=ind2sub([n,n],edges(e))
%     if ii==jj; L(ii,ii)=1; continue; end
%     L(ii,jj)=-1/sqrt(deg(ii)*deg(jj));
% end