function ladle = ladle(X, maxD,nboot)

%Ladle Estimator
%
% INPUT
% X: Data matrix
% maxD: pre-specified maximum number of factors
% nboot: number of bootstrap
%
%OUTPUT
% ladle: selected number of factors by ladle estimator



[n,p] = size(X);
cv = zeros(1, maxD);
np = (p+n)/(p*n);

[W,D,H] = svd(X,0);        
eigenvalue=diag(D);
lambda=eigenvalue(1:maxD+1);
phin=lambda/(1+sum(lambda));


fn0=zeros(maxD+1,1);
for j=1:nboot
   u=datasample(1:n,n);
   Xboot=X(u,:);
   [Wboot,Dboot,Hboot] = svd(Xboot);    
   for r=1:maxD
   fn0(r+1)=fn0(r+1)+1-abs(det(H(:,1:r)'*Hboot(:,1:r)));
   end
end
fn0=fn0/nboot;
fn=fn0/(1+sum(fn0));
gn=phin+fn;

ladle = find(gn == min(gn))-1;  

