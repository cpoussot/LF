%function [Vproj,Wproj,Vproj_x0] = projectors(H,Hr,J,X,Y,la_,mu_,Rj,Li,PROJ)
function [Vproj,Wproj,Vproj_x0] = projectors(H,info_loep,PROJ)

%
Hr  = info_loep.Hr;
J   = info_loep.J;
X   = info_loep.X;
Y   = info_loep.Y;
la_ = info_loep.la;
mu_ = info_loep.mu;
Rj  = info_loep.R;
Li  = info_loep.L;
%
r = length(Hr.a);
n = length(H.a);
if isempty(H.e)
    H.e = eye(n);
end
if isempty(Hr.e)
    Hr.e = eye(r);
end
R = []; O = [];
for ii = 1:numel(la_)
    R = [R  ((la_(ii)*H.e-H.a)\H.b)*Rj(:,ii)];
    O = [O; Li(ii,:)*(H.c/((mu_(ii)*H.e-H.a)))];
end
R   = real(R*J);
O   = real(J'*O);
%norm(-O*H.e*R - LL)
%norm(-O*H.a*R - SS)

% PROJ = 'svd';       % <= original Loewner projected via SVD
% PROJ = 'balreal';   % <= add a balanced realisation projection (for fun)
switch lower(PROJ)
    case 'balreal'
        [Hr,G,T,Ti] = balreal(Hr);
        Vproj       = R*X*Ti;
        Wproj       = (T*Y'*O)';
    case 'svd'
        Vproj       = R*X;
        Wproj       = (Y'*O)';
    case 'none'
        Vproj       = R;
        Wproj       = (O)';
end
Vproj       = Vproj*info_loep.Vchol;
Vproj_x0    = Hr.e\(Wproj');