clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
% MDSPACK
setenv('MDSHOME','/Users/charles/Documents/MDS')
addpath('/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/v1.1.0/API/matlab/')
addpath('/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/v1.1.0/bin')

%%% General variables 
lw      = 3;  % linewidth
mw      = 20; % markersize
CAS     = 2;

%%% Select example
nip = 100;
switch CAS
    case 1 % Random
        ny  = 4;
        nu  = 3;
        S   = rss(20,ny,nu);
        S.E = eye(length(S.A));
        S.D = 2*ones(ny,nu);
    case 2 % Passive RLC
        A   = [-20 -10 0 0 0; 10 0 -10 0 0; 0 10 0 -10 0; 0 0 10 0 -10; 0 0 0 10 -2];
        B   = [20 0 0 0 0]';
        nu  = size(B,2);
        ny  = nu;
        C   = B';
        D   = 1*eye(nu);
        E   = eye(size(A));
        S   = dss(A,B,C,D,E);
    case 3 % Very simple
        ny  = 1;
        nu  = 1;
        S   = ss(tf(1,[1 0.01 1]));
        S.E = eye(length(S.A));
        nip = 4;
    case 4 % Full I/O 
        n   = 7;
        nu  = n;
        ny  = n;
        S   = rss(n,ny,nu);
        S.E = eye(n);
        S.B = eye(n);
        S.C = eye(n);
        S.D = ones(n,n);
end
%
eigS        = eig(S);
[A,B,C,D,E] = dssdata(S);
G           = @(s) C*((s*E-A)\B)+D;

%%% Data
la_ = (logspace(-2,2,nip))*1i;    la_ = sort([la_ conj(la_)]);
mu_ = (logspace(-2,2,nip)+.1)*1i; mu_ = sort([mu_ conj(mu_)]);
la_ = (logspace(-2,2,nip));
mu_ = (logspace(-2,2,nip)+.1);
k   = length(la_);
q   = length(mu_);
R   = 3*ones(nu,k);
L   = ones(q,ny);
for ii = 1:k
    W(1:ny,1:nu,ii) = G(la_(ii));
end
for ii = 1:q
    V(1:ny,1:nu,ii) = G(mu_(ii));
end

%%% Loewner 
opt         = [];
opt.target  = 1e-12;
opt.D       = 0*ones(ny,nu);
opt.real    = false;
% Tangential 
[htng,itng] = lf.loewner_tng(la_,mu_,W,V,R,L,opt);
itng
if ~itng.isCC
    test1 = itng.MU*itng.LL - itng.LL*itng.LA;
    test2 = itng.V*itng.R - itng.L*itng.W;
    norm(test1-test2)
    test1 = itng.MU*itng.SS - itng.SS*itng.LA;
    test2 = itng.MU*itng.V*itng.R - itng.L*itng.W*itng.LA;
    norm(test1-test2)
end
% Block
[hblk,iblk] = lf.loewner_blk(la_,mu_,W,V,opt);
iblk
% Some informations
isreal(itng.Hr)
size(itng.Hr)
isreal(iblk.Hr)
size(iblk.Hr)

%%%
figure, hold on, grid on
plot(itng.sv,'-o','MarkerSize',mw,'LineWidth',lw)
plot(itng.sv_nu,'-x','MarkerSize',mw,'LineWidth',lw)
plot(iblk.sv,'--o','MarkerSize',mw,'LineWidth',lw)
plot(iblk.sv_nu,'--x','MarkerSize',mw,'LineWidth',lw)
set(gca,'YScale','log')
xlabel('$k$','Interpreter','latex')
ylabel('Normalized singular value','Interpreter','latex')
legend({'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)', ... 
        'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)'},'interpreter','latex')

eigHt = eig(itng.Hr);
eigHb = eig(iblk.Hr);
figure, hold on, grid on
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHt),imag(eigHt),'o','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHb),imag(eigHb),'s','MarkerSize',mw,'LineWidth',lw);
xlabel('Real','Interpreter','latex')
ylabel('Imag.','Interpreter','latex')
legend({'Original' 'Loewner tangent' 'Loewner block'},'Interpreter','latex','location','northwest')

w = logspace(-2,3,300);
figure, hold on
mdspack.bodemag(G,w,'-')
mdspack.bodemag(htng,w,'--')
mdspack.bodemag(hblk,w,':')
legend({'Original' 'Loewner tangent' 'Loewner block'},'Interpreter','latex','location','best')

%hdiag  = dss(itng.SSt,itng.Vt,itng.Wt,0,itng.LLt);
% for i = 1:length(itng.lat)
%     Glat(i) = G(itng.lat(i));
%     Hlat(i) = hr(itng.lat(i));
% end
% vpa([Glat(:) Hlat(:) itng.wt(:)],4)
% for i = 1:length(itng.mut)
%     Gmut(i) = G(itng.mut(i));
%     Hmut(i) = hr(itng.mut(i));
% end
% vpa([Gmut(:) Hmut(:) itng.vt(:)],4)
