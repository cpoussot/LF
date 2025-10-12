% Loewner algorithm (tangential passive version)
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [hr,info] = lf.loewner_passive(la,mu,W,V,R,L,D,opt)
%  
% Input arguments
%  - la : interpolation points (k x 1, complex)
%  - mu : interpolation points (q x 1, complex)
%  - W  : k-dimensional ny x nu structure of the function/data 
%         evaluated at points "la" (ny x nu x k, complex)
%         W = H(la), where H(.) is the underlying model
%  - V  : q-dimensional ny x nu structure of the function/data  
%         evaluated at points "mu" (ny x nu x q, complex)
%         V = H(mu), where H(.) is the underlying model
%  - R  : right tangential directions (ny x k, complex)
%  - L  : left tangential directions (q x nu, complex)
%  - D  : H(infty) value (ny x nu, complex)
%  - opt: optional arguments
%    * target: rational order (if integer >=1), 
%              SVD tolerance (if real <1)
%              automatic (if Inf or [])
%    * D     : D-term (ny x nu, complex)
%              /!\ use with discretion
%    * real  : boolean to try (if possible) real realization 
% 
% Output arguments
%  - hr   : approximation model (handle function)
%  - info : structure with informations about the Loewner world
%    * r    : rational order (integer)
%    * nu   : McMilan degree (integer)
%    * isCC : is complex conjugated (boolean)
%             if true, most of the matrices are converted to real
%    * J    : complex to real transformation matrix (k x q, complex)
%    * la   : lambda's column interpolation points (k x 1, complex)
%    * mu   : mu's row interpolation points (q x 1, complex)
%    * sv   : normalized singular values of [LL SS] (min(q,k), real)
%    * sv_nu: normalized singular values of LL (min(q,k), real)
%    * LL   : Loewner matrix (q x k, complex)
%    * SS   : shifted Loewner matrix (q x k, complex)
%    * LA   : Lambda matrix (k x k, complex)
%    * MU   : Mu matrix (q x q, complex)
%    * L,R  : left, right tangential data (data x tangent directions)
%    * V,W  : tangential input and output matrices (q x nu & ny x k, complex)
%    * D    : D-term (ny x nu, complex)
%    * H    : full Loewner form (handle function)
%               H(s)=W(-s LL+SS)\V+D
%    * X,Y  : left and right projectors (k x r & q x r, complex)
%    * Rr,Lr: right and left tangential directions  (nu x r & r x ny, complex)
%    * Hr   : compressed Loewner form (state-space, complex)
%               Hr(s)=Cr(sEr-Ar)\Br+Dr, 
%             where (Er,Ar,Br,Cr,Dr) are available in info.Er ...
%    * lat  : compressed column (right) interpolation points (r x 1, complex)
%    * mut  : compressed row (left) interpolation points (r x 1, complex)
%    * ... Others not documented yet
%    * sz_R     : spectral zeros directions
%    * sz_R_pos : spectral zeros directions positive
%    * sz_la    : spectral zeros 
%    * sz_la_pos: spectral zeros positive
%    * Vchol    : Cholesky decomposition of LL
%    * chol_flag: Cholesky decomposition flag
%    * Hrn      : normalized passive model
% 
% Note 
% Sylvester equations 
%    MU*LL-LL*LA = V*R-L*W 
%    MU*SS-SS*LA = MU*V*R-L*W*LA
% may be checked for complex form only, i.e.
% if ~info.isCC
%   test1 = info.MU*info.LL - info.LL*info.LA;
%   test2 = info.V*info.R - info.L*info.W;
%   norm(test1-test2) % small
%   test1 = info.MU*info.SS - info.SS*info.LA;
%   test2 = info.MU*info.V*info.R - info.L*info.W*info.LA;
%   norm(test1-test2) % small
% end
% 
% Description
% Loewner rules.
%


function [hrp,info] = loewner_passive(la_,mu_,W_,V_,R,L,D,opt)

tol_sz  = 1e-10;
nu      = size(W_,1);
ny      = nu;
%
if nargin < 8 || ~isa(opt,'struct')
    Ds  = 0;
elseif isa(opt,'struct')
    if isfield(opt,'Ds')
        Ds  = opt.Ds*eye(nu);
    else
        Ds  = 0*eye(nu);
    end
end
D   = D  + Ds;
W_  = W_ + Ds;
V_  = V_ + Ds;
% step 1: Loewner classic (with D)
[h_loe,info_loe]    = lf.loewner_tng(la_,mu_,W_-D,V_-D,R,L,opt);
info_loe.Dr         = D;
info_loe.Hr.D       = D;
[info_loe.Hr,~]     = stabsep(info_loe.Hr);
h_loe               = @(s) info_loe.Hr.C*((s*info_loe.Hr.E-info_loe.Hr.A)\info_loe.Hr.B) + info_loe.Hr.D;

% step 2: spectral zeros
[sz_R,sz_la]        = lf.spectral_zeros(info_loe.Hr);
sz_la_pos           = sz_la(real(sz_la)>tol_sz);
sz_R_pos            = sz_R(:,real(sz_la)>tol_sz);

% step 3: spectral zeros Loewner interpolation
clear la R W
for ii = 1:numel(sz_la_pos)
    la(ii)          = sz_la_pos(ii);
    R(:,ii)         = sz_R_pos(end-nu+1:end,ii)/norm(sz_R_pos(end-nu+1:end,ii));
    W(1:ny,1:nu,ii) = h_loe(sz_la_pos(ii));
end
%
%opt.stable      = true;
opt.D           = D;
opt.target      = -1;
[~,info2]       = lf.loewner_tng(la,-conj(la),W,-conj(W),R,R',opt);
info2.Hr.D      = info2.Hr.D - Ds;
%tol_hermite = 1e-10;
%[isstable(info2.Hr) isPassive(info2.Hr) norm(info2.LL-conj(info2.LL.'))<tol_hermite]
% step 4: normalized passive
[T,flag]    = chol(info2.LL); % R'*R = E (E > 0)
if flag == 0
    invT    = T\eye(length(T));
else
    warning('Cholesky decomposition issue')
    invT    = eye(length(info2.Hr.A));
end
Hr          = ss(-invT'*info2.Hr.A*invT,-invT'*info2.Hr.B,info2.Hr.C*invT,info2.Hr.D);
rr          = length(Hr.A);
In          = eye(rr);
hrp         = @(s) Hr.C*((s*In-Hr.A)\Hr.B)+Hr.D;
%
info            = info2;
info.X          = eye(rr);
info.Y          = eye(rr);
% 
info.sz_R       = sz_R;
info.sz_R_pos   = sz_R_pos;
info.sz_la      = sz_la;
info.sz_la_pos  = sz_la_pos;
%
info.Vchol      = invT;
info.chol_flag  = flag;
info.Hrn        = Hr;
