% Based on: Sales-Pardo et al, "Extracting the hierarchical organization 
%                of complex systems", PNAS, Sep 25, 2007; vol.104; no.39 
% Supplementary material: 
% http://www.pnas.org/content/suppl/2008/02/27/0703740104.DC1/07-03740SItext.pdf
%
% INPUTs: N: number of nodes; L: number of hierarchy levels; 
%         [G1,G2,..,GL]: number of nodes in each group in each level
%         kbar: average degree, 
%         rho [optional]: ratio between average degrees at different levels
%                                                   (see supplementary material) 
%         Example inputs (from paper): N=640, L=3, G=[10,40,160], kbar=16, rho=1
% OUTPUTs: edge list, in mx2 or mx3 format, where m = number of edges
%
% Other routines used: symmetrizeEdgeL.m
% GB: last updated, November 24 2012


function eL = nestedHierarchiesModel(N,L,G,kbar,rho)

% first check whether the inputs are of the right size/type

if length(G)~=L; printf('The number of levels do not match. length(G) should be L');  return; end

for x=2:L
    if G(x)/G(x-1)~=ceil(G(x)/G(x-1)); printf('number of groups not an integer at level %2i\n',x); return ; end
end
if N/G(L)~=ceil(N/G(L)); printf('number of groups not an integer at level %2i\n',L); return ; end


% formula on page 3 of supplementary material
if nargin<5; rho = kbar/(G(L)-1) - 1 + 0.05; end;  % set to lower bound

if rho < kbar/(G(L)-1) - 1; printf('rho is below its theoretical lower bound, given kbar and G(L)\n'); return; end
% ..........................................................


% create node membership to various nested groups 

belongsto = {};
for ii=1:N; belongsto{ii} = zeros(1,L);  end % the level 1, 2, 3,... groups ii belongs to


for ii=1:N  % across all nodes
  
  for level = 1:L   % across all levels
    
    group = (ii-mod(ii,G(level)))/G(level) + 1;
    if mod(ii,G(level))==0; group -= 1; end
    
    belongsto{ii}(level) = group;
    
  end
end
% ..........................................................


eL = [];

for i=1:N
    for j=i+1:N
        
        x = sum( belongsto{i}==belongsto{j} );   % number of common groups i and j belong to
        
        if x==0
            pij = (rho/(1+rho))^(L-x) * kbar / (N - G(1));            
        else
            pij = rho^(L-x)/(1+rho)^(L-x+1) * kbar/(G(x)-1);        
        end
        
        if rand < pij; eL = [eL; i j 1]; end
        
        
    end
end

eL = symmetrizeEdgeL(eL);