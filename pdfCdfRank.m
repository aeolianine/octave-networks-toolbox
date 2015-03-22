% Compute the pdf, cdf and rank distributions for a sequence of numbers.
%
% INPUTS: sequence of values: x, (1xn), 'plt' - 'on' or 'off'
%         bin - bin size, default is [], then it gets selected
%         automatically
% OUTPUTS: pdf, cdf and rank distribution values, plot is optional
%
% Note: pdf = frequency distribution, cdf = cumulative frequency, 
%       rank = log-log scale of the sorted sequence
% GB: last updated, November 24 2012

function [xpdf,ypdf,xcdf,ycdf,logk,logx]=pdfCdfRank(x,plt='off',bin=[])

xx=unique(x);

if length(bin)==0; bin = length(xx)/5; end


for ii=1:numel(xx)
  xcdf(ii) = xx(ii);
  ycdf(ii) = length(find(x<=xx(ii)))/numel(x);
  
  % how many x's fall in the interval [xx(ii)-0.5*numel(xx)/bin,xx(ii)+0.5*numel(xx)/bin]
  xpdf(ii) = xx(ii);
  ypdf(ii) = length(find(abs(xx(ii)-x)<=0.5*numel(xx)/bin))/numel(x); 
end

x=-sort(-x);
logk=log(1:length(x));
logx=log(x);

if strcmp(plt,'on')
     set(gcf,'color',[1,1,1])
     
     subplot(2,2,1)
     hist(x,bin,'facecolor',[0.5 0.5 0.5])
     title('histogram')
     
     subplot(2,2,3)
     plot(xpdf,ypdf,'k.')
     title('pdf')
     axis('tight')
     
     subplot(2,2,2)
     plot(xcdf,ycdf,'k.')
     title('cdf')
     axis('tight')
     
     subplot(2,2,4)
     plot(logk,logx,'k.')
     title('rank')
     axis('tight')
end