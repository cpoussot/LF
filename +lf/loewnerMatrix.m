function [LL,SS] = loewnerMatrix(la_,mu_,W,V)

[ny,nu] = size(W(:,:,1));
ni      = length(mu_);
nj      = length(la_);
LL      = zeros(ni*ny,nj*nu);
SS      = zeros(ni*ny,nj*nu);
co0     = 1;
for jj = 1:nj
    li0 = 1;
    co  = co0:co0+nu-1; 
    for ii = 1:ni
        li        = li0:li0+ny-1;
        LL(li,co) = (reshape(V(:,:,ii),ny,nu)-reshape(W(:,:,jj),ny,nu))/(mu_(ii)-la_(jj));
        SS(li,co) = (mu_(ii)*reshape(V(:,:,ii),ny,nu)-reshape(W(:,:,jj),ny,nu)*la_(jj))/(mu_(ii)-la_(jj));
        li0       = li0+ny;
    end
    co0 = co0+nu;
end
