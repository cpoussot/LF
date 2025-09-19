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
lw  = 4;  % linewidth
mw  = 20; % markersize

%%% Select example
% 'siso_passive_simple'
% 'siso_passive_aca'
% 'siso_passive_gugercin'
[G,S,nu,ny]     = lf.examples('siso_passive_gugercin');
[SZ_v,SZ,~,~]   = lf.spectral_zeros(S);
r = 2%1e-12

%%% Data
% nip     = 50;
% puls    = logspace(-1,2,nip);
% la      = puls(1:2:end)*1i; la = sort([la conj(la)]);
% mu      = puls(2:2:end)*1i; mu = sort([mu conj(mu)]);
pts     = linspace(.1,50,20);
la      = pts(1:2:end);
mu      = pts(2:2:end);
k       = length(la);
q       = length(mu);
R       = ones(nu,k);
L       = ones(q,ny);
for ii = 1:k
    W(1:ny,1:nu,ii) = G(la(ii));
end
for ii = 1:q
    V(1:ny,1:nu,ii) = G(mu(ii));
end

%%% Loewner 
opt         = [];
opt.target  = r;
opt.real    = true;
[htng,itng] = lf.loewner_tng(la,mu,W,V,R,L,opt);
isPassive(itng.Hr)
%[SZr_v,SZr,~,~] = lf.spectral_zeros(itng.Hr);

%%% Loewner pH
% step 1: Loewner
D           = S.D;
opt         = [];
opt.target  = r;
opt.real    = true;
[h1,i1]     = lf.loewner_tng(la,mu,W-D,V-D,R,L,opt);
i1.Dr       = D;
i1.Hr.D     = D;
h1          = @(s) i1.Hr.C*((s*i1.Hr.E-i1.Hr.A)\i1.Hr.B) + i1.Hr.D;
[isstable(i1.Hr) isPassive(i1.Hr)]
% step 2: spectral zeros
m = nu;
[sz_R,sz_la]    = lf.spectral_zeros(i1.Hr);
tol             = 1e-10;
sz_la_pos       = sz_la(real(sz_la)>tol);
sz_R_pos        = sz_R(:,real(sz_la)>tol);
%
figure, hold on, grid on
plot(real(SZ),imag(SZ),'.','MarkerSize',mw,'LineWidth',lw)
plot(real(sz_la),imag(sz_la),'o','MarkerSize',mw,'LineWidth',lw)
plot(real(sz_la_pos),imag(sz_la_pos),'s','MarkerSize',mw,'LineWidth',lw)
title(['Spectral zeros $r=' num2str(r) '$'])
legend({'Original' 'Loewner' 'Loewner (positive only)' 'passive Loewner'},'location','best')

% step 3: Loewner interpolate positive spectral zeros
clear R W la
for ii = 1:numel(sz_la_pos)
    la(ii)          = sz_la_pos(ii);
    R(:,ii)         = sz_R_pos(end-m+1:end,ii)/norm(sz_R_pos(end-m+1:end,ii));
    W(1:ny,1:nu,ii) = h1(sz_la_pos(ii));
end
%
opt.D       = D;
opt.real    = true;
opt.target  = -1;
[hp,ip]     = lf.loewner_tng(la,-conj(la),W,-conj(W),R,conj(R).',opt);
[isstable(ip.Hr) isPassive(ip.Hr)]
ishermitian(ip.LL)

%%%
figure, hold on, grid on
plot(itng.sv,'-o','MarkerSize',mw,'LineWidth',lw)
plot(itng.sv_nu,'-o','MarkerSize',mw,'LineWidth',lw)
plot(ip.sv,'--x','MarkerSize',mw,'LineWidth',lw)
plot(ip.sv_nu,'--x','MarkerSize',mw,'LineWidth',lw)
set(gca,'YScale','log')
xlabel('$k$'), ylabel('Normalized singular value')
title(['Singular values $r=' num2str(r) '$'])
legend({'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)' ...
        'svd($[\bf{L},\bf{M}]$) pH','svd($\bf{L}$) pH'},'Location','best')

eigS    = eig(S);
eigHt   = eig(itng.Hr);
eigHph  = eig(ip.Hr);
figure, hold on, grid on
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHt),imag(eigHt),'o','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHph),imag(eigHph),'s','MarkerSize',mw,'LineWidth',lw);
title(['Eigenvalues $r=' num2str(r) '$'])
xlabel('Real'), ylabel('Imag.')
legend({'Original' ['Loewner (\texttt{isPassive}:' num2str(isPassive(itng.Hr)) ')'] ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(ip.Hr)) ')']},'location','best')
%
if isPassive(ip.Hr)
    %%% pH form in real struture
    [A,B,C,D]       = ssdata(ip.Hr);
    Sblk            = [-A -B; C D];
    n               = length(A);
    JGN             = 1/2*(Sblk-Sblk');
    RPS             = 1/2*(Sblk+Sblk');
    info.J          = -JGN(1:n,1:n);
    info.G          = -JGN(1:n,n+1:end);
    info.N          = JGN(n+1:end,n+1:end);
    info.R          = RPS(1:n,1:n);
    info.P          = RPS(1:n,n+1:end);
    info.S          = RPS(n+1:end,n+1:end);
    info.Q          = eye(n);
    In              = eye(n);
    hph             = @(s) ((info.G+info.P)')*((s*In-(info.J-info.R))\(info.G-info.P))+(info.N+info.S);
    info.chol_flag  = flag;
    Wcal            = [info.R info.P; info.P' info.S];
    Vcal            = [-info.J -info.G; info.G' info.N];

    [ishermitian(Wcal) ishermitian(Vcal,'skew') all(eig(info.Q)>0)]
end

w = logspace(-2,3,1e4);
figure, hold on
mdspack.bodemag(G,w,'-')
mdspack.bodemag(htng,w,'--')
mdspack.bodemag(hp,w,':')
mdspack.bodemag(hph,w,'-.')
title(['Bode gain $r=' num2str(r) '$'])
legend({'Original' ['Loewner (\texttt{isPassive}:' num2str(isPassive(itng.Hr)) ')'] ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(ip.Hr)) ')'] 'pH Loewner'},'location','best')
figure, hold on
mdspack.bodephase(G,w,'-')
mdspack.bodephase(htng,w,'--')
mdspack.bodephase(hp,w,':')
mdspack.bodephase(hph,w,'-.')
title(['Bode phase $r=' num2str(r) '$'])
legend({'Original' ['Loewner (\texttt{isPassive}:' num2str(isPassive(itng.Hr)) ')'] ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(ip.Hr)) ')'] 'pH Loewner'},'location','best')

figure, hold on
mdspack.nyquist(G,w,'-')
mdspack.nyquist(htng,w,'--')
mdspack.nyquist(hp,w,':')
mdspack.nyquist(hph,w,'-.')
title(['Nyquist $r=' num2str(r) '$'])
legend({'Original' ['Loewner (\texttt{isPassive}:' num2str(isPassive(itng.Hr)) ')'] ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(ip.Hr)) ')'] 'pH Loewner'},'location','best')