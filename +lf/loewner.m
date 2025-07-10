% Syntax
% [hr,Hr,sv,LL,SS] = loewner(la,mu,W,V,R,L,robj,goreal)
%  
% Input arguments
%  - la     : interpolation points (n1 x 1, complex)
%  - mu     : interpolation points (n2 x 1, complex)
%  - W      : n1-dimensional ny x nu structure of the function/data 
%             evaluated at points "la" (ny x nu x n1, complex)
%             W = H(la), where H(.) is the underlying model
%  - V      : n2-dimensional ny x nu structure of the function/data  
%             evaluated at points "mu" (ny x nu x n2, complex)
%             V = H(mu), where H(.) is the underlying model
%  - R      : right tangential directions (ny x n1, complex)
%  - L      : left tangential directions (n2 x nu, complex)
%  - robj   : rational order (integer > 0)
%  - goreal : enforce realness (boolean)
% 
% Output arguments
%  - hr     : approximation (handle function)
%  - Hr     : approximation (state-space form)
%  - sv     : normalized singular values
%  - LL     : Loewner matrix (n2 x n1 real)
%  - LL     : shifted Loewner matrix (n2 x n1 real)
%  - info   : structure with a lot of informations about the Loewner world
%               r      : rational order
%               la     : lambda's (column interpolation points)
%               mu     : mu's (row interpolation points)
%               sv     : normalized singular values of [LL SS]
%               LL     : Loewner matrix
%               SS     : shifted Loewner matrix
%               lar    : compressed column (right) interpolation points
%               mur    : compressed row (left) interpolation points
%               E,A,B,C: H(s)=C(sE-A)\B
% 
% Description
% Loewner rules
% 

function [H,Hr,TFr,eigHr,LL,SS] = loewner(la_,mu_,W_,V_,R,L,robj)

k           = length(la_);
q           = length(mu_);
[ny,nu,~]   = size(W_);
%
for ii = 1:k
    W(1:ny,ii) = W_(:,:,ii)*R(:,ii);
end
for ii = 1:q
    V(ii,1:nu) = L(ii,:)*V_(:,:,ii);
end
for ii = 1:q
    for jj = 1:k
        num1        = V(ii,:)*R(:,jj) - L(ii,:)*W(:,jj);
        num2        = mu_(ii)*V(ii,:)*R(:,jj) - L(ii,:)*W(:,jj)*la_(jj);
        den         = mu_(ii)-la_(jj);
        LL(ii,jj)   = num1/den;
        SS(ii,jj)   = num2/den;
    end
end
% Go real
J0  = (1/sqrt(2))*[1 1i; 1 -1i];
J   = [];
kk  = 1;
while length(J) < length(la_)
    if imag(la_(kk)) == 0
        J   = blkdiag(J,1);
        kk  = kk + 1;
    else
        J   = blkdiag(J,J0);
        kk  = kk + 2;
    end
end
LL  = real(J'*LL*J);
SS  = real(J'*SS*J);
V   = real(J'*V);
W   = real(W*J);

%%% Truncate
[L1,S1,R1]  = svd([LL,SS],'econ');
[L2,S2,R2]  = svd([SS',LL']','econ');
r           = min(sum(diag(S1)/S1(1,1)>1e-12),robj);
if isempty(robj)
    r   = length(LL);
end
%
Y         	= L1(:,1:r);
X          	= R2(:,1:r);
%
H          	= dss(-SS,V,W,0,-LL);
Hr          = dss(-Y'*SS*X,Y'*V,W*X,0,-Y'*LL*X);
TFr         = @(s) Hr.c*((Hr.e*s-Hr.a)\Hr.b) + Hr.d;
eigHr       = eig(Hr);


% %%% Projection
% On  = ones(q,1);
% Y   = L1(:,1:r);
% X   = R2(:,1:r);
% Ar  = -Y'*SS*X;
% Br  = Y'*V;
% Cr  = W*X;
% Er  = -Y'*LL*X;
% Lr  = Y'*On;
% Rr  = On.'*X;
% %%% Compressed IP
% evR = eig(Er\(Ar+Br*Rr));
% evL = eig((Ar+Lr*Cr)/(Er));
% %%%
% info.r      = length(Ar);
% info.la     = la;
% info.mu     = mu;
% info.sv     = sv;
% info.LL     = LL; 
% info.SS     = SS;
% info.lar    = evR;
% info.mur    = evL;
% %
% info.E      = Er;
% info.A      = Ar;
% info.B      = Br;
% info.C      = Cr;