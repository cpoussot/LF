% Loewner block matrix
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
%
% Syntax
%  [LL, SS] = lf.loewnerMatrix(la, mu, W, V)
%  
% Input arguments
%  - la : interpolation points (k x 1, complex)
%  - mu : interpolation points (q x 1, complex)
%  - W  : k-dimensional ny x nu structure of the function/data evaluated 
%         at points "la" (ny x nu x k, complex)
%         W = H(la), where H(.) is the underlying model
%  - V  : q-dimensional ny x nu structure of the function/data evaluated 
%         at points "mu" (ny x nu x q, complex)
%         V = H(mu), where H(.) is the underlying model
% 
% Output arguments
%  - LL : MIMO Loewner matrix (q.ny x k.nu, complex)
%  - SS : MIMO shifted Loewner matrix (q.ny x k.nu, complex)
% 
% Description
% Computes the MIMO Loewner matrix in the single variable case, using the
% system's data obtained from interpolation points, without using 
% tangential interpolation. The data are as
%  W = H(la) and V = H(mu)
% where H(s) is the single variable function to interpolate.
%

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
