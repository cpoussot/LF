clearvars; close all; clc;

%%% General variables 
lw      = 3;  % linewidth
mw      = 20; % markersize
%
robj    = 1%1e-12;
S       = ss(tf(1,[1 0.2 1]))
nu = 1; ny = nu;

d = 0%[];
eigS        = eig(S);
[A,B,C,D,E] = dssdata(S);
G           = @(s) C*((s*E-A)\B)+D;

%%% Data
la_ = [1 2 3];
mu_ = -la_;
k   = length(la_);
q   = length(mu_);
R   = ones(nu,k);
L   = ones(q,ny);
for ii = 1:k
    W(1:ny,1:nu,ii) = G(la_(ii));
end
for ii = 1:q
    V(1:ny,1:nu,ii) = G(mu_(ii));
end

%%% Loewner pencil
[hr,info]  = lf.loewner(la_,mu_,W,V,R,L,robj,d);
info

for i = 1:length(info.lar)
    Glar(i) = G(info.lar(i));
    Hlar(i) = hr(info.lar(i));
end
abs(Glar(:)-Hlar(:))
for i = 1:length(info.mur)
    Gmur(i) = G(info.mur(i));
    Hmur(i) = hr(info.mur(i));
end
abs(Gmur(:)-Hmur(:))

%sym(info.Er) % LL
%sym(info.Ar) % LLs
%sym(info.Br) % V
%sym(info.Cr) % W

%eig(inv(info.Er)*(info.Ar+info.Br*info.R))

% figure, hold on, grid on
% plot(info.sv,'-o','MarkerSize',mw,'LineWidth',lw)
% plot(info.sv_nu,'-x','MarkerSize',mw,'LineWidth',lw)
% set(gca,'YScale','log')
% xlabel('$k$','Interpreter','latex')
% ylabel('Normalized singular value','Interpreter','latex')
% legend({'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)'},'interpreter','latex')
% 
% eigH = eig(info.Hr);
% figure, hold on, grid on
% plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw);
% plot(real(eigH),imag(eigH),'o','MarkerSize',mw,'LineWidth',lw);
% xlabel('Real','Interpreter','latex')
% ylabel('Imag.','Interpreter','latex')
% legend({'Original' 'Loewner'},'Interpreter','latex','location','northwest')
% 
