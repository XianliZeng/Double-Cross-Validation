function [num_factor, dcv] = DoubleCV(X, K, maxD)

% K-fold DCV
%
% INPUT
% X: Data matrix  
% K: the number of folds, negagive means leave-one-out. e.g. K=5 means
%        5-fold CV; 
% maxD: pre-specified maximum number of factors, it can be any number
%        smaller than p. It is only set for reducing the unnecessary
%        calculations
%
%OUTPUT
% num_factor: selected number of factors by DCV
% dcv:  dcv values for working number of factors: 0,1,...,K

[n,p] = size(X);

if p > n
    X = X';
    p0 = n;
    n = p;
    p = p0;
end

[~, rp] = sort(rand(n, 1));
X=X(rp,:);
if K <=1; K = n; end;
K = min(K, n);

startsEnd = zeros(K, 2);
nK = floor(n/K);
startsEnd(:,2) = nK;
m = n - nK*K;
startsEnd(1:m,2) = nK+1;
startsEnd(:,2) = cumsum(startsEnd(:,2));
startsEnd(:,1) = 1 + [0; startsEnd(1:(K-1),2)];


cv = zeros(K,maxD+1);    
for k = 1:K
    Ik = startsEnd(k,1):startsEnd(k,2);
    Jk = setdiff(1:n, Ik);
    
    Xi = X(Jk,:);
    yi = X(Ik,:)';
    if p > n
        [H,D,W] = svd(Xi', 0); 
    else
        [W,D,H] = svd(Xi, 0);
    end
    zi = H'*yi;
    S = (1-cumsum(H.^2,2)).^2;
    for d = 0:maxD
        if d > 0
            e = sum(sum((yi-H(:,1:d)*zi(1:d,:)).^2, 2)./S(:,d));  
        else
            e = sum(sum((yi-repmat(mean(yi), size(yi,1), 1)).^2));
        end
        cv(k, d+1) = e;
    end
end
dcv = mean(cv)/nK/n;
%cv = log(cv) + (0:maxD)/np*(log(np));
num_factor = find(dcv == min(dcv))-1;    
