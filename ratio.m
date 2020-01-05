function [rationumber] = ratio(X, maxD)


% Ratio-based test
%
% INPUT
% X: Data matrix

% maxD: pre-specified maximum number of factors
%
%OUTPUT
% rationumber: selected number of factors by Ratio-based test

[n,p] = size(X);
ratio=zeros(maxD,1);


[W,D,H] = svd(X);   
D=diag(D);
for d = 1:maxD
   ratio(d)=D(d+1)/D(d);
end

rationumber = find(ratio == min(ratio));  
