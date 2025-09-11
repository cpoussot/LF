function [SZ_v,SZ,AA,EE] = spectral_zeros(H)

A   = H.A;
B   = H.B;
C   = H.C;
D   = H.D;
E   = H.E;
n   = length(A);
m   = size(B,2);
if isempty(E)
    E = eye(n);
end
AA          = [zeros(n,n), A, B; ...
               A', zeros(n,n) C'; ...
               B', C, D + D'];
EE          = [zeros(n,n) E zeros(n,m); ...
               -E', zeros(n,n), zeros(n,m); ...
               zeros(m,n), zeros(m,n),zeros(m,m)];
[SZ_v,SZ]   = eig(AA,EE);
SZ          = diag(SZ);
%SZ(1:end-m) % Select finite eigenvalues
%SZ_v(:,1:end-m) % Select finite eigenvalues

% Select finite eigenvalues
n_discard = 0;
if isnan(norm(SZ)) || isinf(norm(SZ)) || any(abs(SZ) > 1e10)
    idxNaN          = find(isnan(SZ));
    idxInf          = find(isinf(SZ));
    idx1e10         = find(abs(SZ)>1e10);
    discard         = unique([idxNaN; idxInf; idx1e10]);
    n_discard       = numel(discard);
    keep            = setdiff(1:numel(SZ),discard);
    SZ              = SZ(keep);
    %r               = SZ_v(:,discard);
    SZ_v            = SZ_v(:,keep);
    %fprintf('>> Infinite/NaN/>1e10 spectral zeros removed: %d\n',n_discard)
end
fprintf('>> Sectral Zeros Finite   : %d (th), %d (num) \n',2*n,numel(SZ))
fprintf('>>               Infinite : %d (th), %d (num) \n',m,n_discard)

% % Check spectral zeros solves the Z(s_j)r_j + Z^T(-s_j)r_j = 0
% Hs  = @(s) C*((s*E-A)\B)+D;
% for jj = 1:numel(SZ)
%     sj      = SZ(jj);
%     rj      = r(jj,:).';
%     err(jj) = norm(Hs(sj)*rj + (Hs(-sj)).'*rj);
% end
% if max(err) > 1e-6
%     check = sum(err>1e-6);
%     fprintf('>> Spectral condition Z(s_j)r_j + Z^T(-s_j)r_j = 0 not validated, %d times\n',check)
% end    