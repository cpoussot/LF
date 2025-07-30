% Loewner algorithm (block version)
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
%
% Syntax
% [hr,info] = lf.loewner_blk(la,mu,W,V,opt)
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
%  - opt: optional arguments
%    * target: rational order (if integer >=1), 
%              SVD tolerance (if real <1)
%              automatic (if Inf or [])
%    * D     : D-term (ny x nu, complex)
%              /!\ use with discretion
% 
% Output arguments
%  - hr   : approximation model (handle function)
%  - info : structure with informations about the Loewner world
%    * r    : rational order (integer)
%    * nu   : McMilan degree (integer)
%    * isCC : is complex conjugated (boolean)
%             if true, most of the matrices are converted to real
%    * la   : lambda's column interpolation points (k x 1, complex)
%    * mu   : mu's row interpolation points (q x 1, complex)
%    * sv   : normalized singular values of [LL SS] (min(ny.q,nu.k), real)
%    * sv_nu: normalized singular values of LL (min(ny.q,nu.k), real)
%    * LL   : Loewner matrix (ny.q x nu.k, complex)
%    * SS   : shifted Loewner matrix (ny.q x nu.k, complex)
%    * V,W  : tangential input and output matrices (ny.q x nu & ny x nu.k, complex)
%    * D    : D-term (ny x nu, complex)
%    * H    : full Loewner form (handle function)
%               H(s)=W(-s SS+LL)\V+D
%    * X,Y  : left and right projectors (nu.k x r & ny.q x r, complex)
%    * Hr   : compressed Loewner form (state-space, complex)
%               Hr(s)=Cr(sEr-Ar)\Br+Dr, 
%             where (Er,Ar,Br,Cr,Dr) are available in info.Er ...
%    * ... Others not documented yet
% 
% Description
% Loewner rules.
%

function [hr,info] = loewner_blk(la_,mu_,W,V,opt)

TOL_SV  = 1e-13;
TOL_CC  = 1e-13;
%
[ny,nu,~]   = size(W);
if nargin < 5 || ~isa(opt,'struct')
    D    = zeros(ny,nu);
    robj = inf;
elseif isa(opt,'struct')
    if isfield(opt,'target')
        robj = opt.target;
    else
        robj = inf;
    end
    if isfield(opt,'D')
        D = opt.D;
    else
        D = zeros(ny,nu);
    end
end
%
k   = length(la_);
q   = length(mu_);
% 
isCC = false;
if (abs(sum(imag(la_.')))<TOL_CC) && ...
   (abs(sum(imag(mu_.')))<TOL_CC) && ...
   (q==k)
    isCC = true;
end

%%% Reshape data
VV  = zeros(q*ny,nu);
li  = 1;
for ii = 1:q
    li       = li:li+ny-1;
    VV(li,:) = reshape(V(:,:,ii),[size(V,1) size(V,2)]);
    li       = li+ny;
end
WW  = zeros(ny,k*nu);
co  = 1;
for ii = 1:k
    co       = co:co+nu-1;
    WW(:,co) = reshape(W(:,:,ii),[size(W,1) size(W,2)]);
    co       = co+nu;
end

%%% Loewner matrices
[LL,SS] = lf.loewnerMatrix(la_,mu_,W,V);

%%% D-term
if ~norm(D) == 0
    [qq,kk] = size(LL);
    L       = ones(qq,ny);
    R       = ones(nu,kk);
    SS      = (SS - L*D*R);
    VV      = VV - L*D;
    WW      = WW - D*R;
end

%%% Go real
if isCC
    J0  = 1/sqrt(2) * [1 1i; 1 -1i]; % J0'*J0 = I
    Jl  = kron(eye(length(la_)/2),kron(J0,eye(ny)));
    Jr  = kron(eye(length(mu_)/2),kron(J0,eye(nu)));
    LL  = real(Jl'*LL*Jr);
    SS  = real(Jl'*SS*Jr);
    WW  = real(WW*Jr);
    VV  = real(Jl'*VV);
end

%%% Compressed model
% orders
[L1,S1,~]   = svd([LL,SS],'econ','vector');
[~,S2,R2]   = svd([SS',LL']','econ','vector');
S_nu        = svd(LL,'econ','vector');
sv          = S1/S1(1,1);
sv_nu       = S_nu/S_nu(1,1);
nu_         = sum(sv_nu>TOL_SV);
if isempty(robj) | isinf(robj)
    r   = sum(sv>TOL_SV);
elseif robj < 1
    r   = sum(sv>robj);
elseif robj >= 1
    r   = robj;
end
Y   = L1(:,1:r);
X   = R2(:,1:r);
% compressed
Er  = -Y'*LL*X;
Ar  = -Y'*SS*X;
Br  = Y'*VV;
Cr  = WW*X;
Dr  = D;
%Lr  = Y'*L;
%Rr  = R*X;

%%% Output
hr  = @(s) Cr*((Er*s-Ar)\Br) + Dr;
h   = @(s) W*((-LL*s+SS)\V) + D;
Hr  = dss(Ar,Br,Cr,Dr,Er);

% %%% Compressed IP
% [Tla,LAt,~] = eig(Ar+Br*Rr,Er); % right IP
% [~,MUt,Tmu] = eig(Ar+Lr*Cr,Er); % left IP
% %[Tla,LAt] = eig(inv(Er)*(Ar+Br*Rr));
% %[Tmu,MUt] = eig((Ar+Lr*Cr)*inv(Er));
% LAt         = diag(LAt);
% MUt         = diag(MUt);
% LLt         = Tmu*Er*Tla;
% SSt         = Tmu*Ar*Tla;
% Vt          = Tmu*Br; % "-"?
% Wt          = Cr*Tla; % "-"?
% Lt          = Tmu*Lr;
% Rt          = Rr*Tla;

%%% Information
% Loewner
info.r      = r;
info.nu     = nu_;
info.isCC   = isCC;
info.la     = la_(:);
info.mu     = mu_(:);
info.sv     = sv;
info.sv_nu  = sv_nu;
info.LL     = LL; 
info.SS     = SS;
info.V      = VV;
info.W      = WW;
info.D      = D;
info.H      = h;
% Compression
info.X      = X;
info.Y      = Y;
% info.Rr     = Rr;
% info.Lr     = Lr;
info.Hr     = Hr;
info.Er     = Er;
info.Ar     = Ar;
info.Br     = Br;
info.Cr     = Cr;
info.Dr     = Dr;
% % Barycentric/modal form
% info.lat    = LAt;
% info.mut    = MUt;
% if (ny == 1) && (nu == 1) 
%     info.vt = Vt./Lt; %-(Tmu*Y'*V)./(Tmu*Y'*L);
%     info.wt = Wt./Rt; %-(W*X*Tla)./(R*X*Tla);
% end
% info.LLt    = LLt; 
% info.SSt    = SSt;
% info.Vt     = Vt;
% info.Wt     = Wt;
