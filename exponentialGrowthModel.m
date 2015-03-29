% Construct a graph with an exponential degree distribution.
% Probability of node s having k links at time t: 
%    p(k,s,t)=1/t*p(k-1,s,t-1)+(1-1/t)*p(k,s,t-1)
%
% INPUTS: number of time steps, t
% OUTPUTs: edge list, mx3
%
% GB, last updated: Nov 11, 2012

function el=exponentialGrowthModel(t)

el=[1 2 1; 2 1 1]; % initialize with two connected nodes

% for all remaining time t
for i=3:t; r = randi(i-1); el=[el; i r 1; r i 1]; end