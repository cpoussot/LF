clearvars; close all; clc;

%%% General variables 
lw      = 3;  % linewidth
mw      = 20; % markersize
CAS     = 1;
%
robj    = inf;

%%% Select example
switch CAS
    case 0 % Random
        ny  = 1;
        nu  = 1;
        S   = rss(20,ny,nu);
        S   = S+1.01*norm(S,inf);
        S.E = eye(length(S.A));
    case 1 % Passive RLC
        A   = [-20 -10 0 0 0; 10 0 -10 0 0; 0 10 0 -10 0; 0 0 10 0 -10; 0 0 0 10 -2];
        B   = [20 0 0 0 0]';
        nu  = size(B,2);
        ny  = nu;
        C   = B';
        D   = eye(nu);
        E   = eye(size(A));
        S   = dss(A,B,C,D,E);
end
eigS        = eig(S);
[A,B,C,D,E] = dssdata(S);
G           = @(s) C*((s*E-A)\B)+D;

%%% Data
la_ = (logspace(0,2,50))*1i;      la_ = sort([la_ conj(la_)]);
mu_ = (logspace(0,2,50)+.1)*1i;   mu_ = sort([mu_ conj(mu_)]);
k   = length(la_);
q   = length(mu_);
R   = ones(nu,k);
L   = ones(q,ny);
for ii = 1:k
    W_(1:ny,1:nu,ii) = G(la_(ii));
end
for ii = 1:q
    V_(1:ny,1:nu,ii) = G(mu_(ii));
end

%%% Loewner pencil
[H,Hr,TFr,eigHr,LL,SS]  = lf.loewner(la_,mu_,W_,V_,R,L,robj);
size(S)
size(H)
size(Hr)
%
figure, hold on, grid on
sv    = svd([LL,SS]);
sv_nu = svd(LL);
plot(sv/sv(1),'-o','MarkerSize',mw,'LineWidth',lw)
plot(sv_nu/sv_nu(1),'-x','MarkerSize',mw,'LineWidth',lw)
set(gca,'YScale','log')
xlabel('$k$','Interpreter','latex')
ylabel('Normalized singular value','Interpreter','latex')
legend({'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)'},'interpreter','latex')

figure, hold on, grid on
plot(eigS,'.','MarkerSize',mw,'LineWidth',lw);
plot(eigHr,'o','MarkerSize',mw,'LineWidth',lw);
xlabel('Real','Interpreter','latex')
ylabel('Imag.','Interpreter','latex')
legend({'Original' 'Loewner'},'Interpreter','latex','location','northwest')

figure, hold on
bodemag(S,'-',Hr,'r--',{1e-2,1e3})
legend({'Original' 'Loewner'},'Interpreter','latex','location','northwest')

