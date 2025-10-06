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
