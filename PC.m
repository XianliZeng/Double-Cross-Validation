function [DPC1,DIC1] = PCp1(X, maxD)


% Penal Criterion
%
% INPUT
% X: Data matrix

% maxD: pre-specified maximum number of factors
%
%OUTPUT
% DPC1: selected number of factors by PC_1
% DIC1: selected number of factors by IC_1



[n,p] = size(X);
cv = zeros(1, maxD);
np = (p+n)/(p*n);

[W,D,H] = svd(X);        
H = H';
Wmax = W(:,1:maxD);
Dmax = D(1:maxD,1:maxD);
Hmax = H(1:maxD,:); 
sigmahat=mean(mean((X - Wmax*Dmax*Hmax).^2));
for d = 0:maxD
    if d > 0
         Wi = W(:,1:d);
         Di = D(1:d,1:d);
         Hi = H(1:d,:);        
        s2 = mean(mean((X - Wi*Di*Hi).^2));
        PC1(d+1)=s2+sigmahat*d*np*log(1/np);

        IC1(d+1) = log(s2) + d*np*log(1/np);
    
    else
        s2 = mean(mean((X - repmat(mean(X, 2), 1, p)).^2));
        PC1(d+1)=s2;

        IC1(d+1) = log(s2);

    end
end

DPC1 = find(PC1 == min(PC1))-1;  
DIC1 = find(IC1 == min(IC1))-1;

