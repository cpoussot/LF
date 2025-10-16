clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end

%%% General variables 
lw  = 4;  % linewidth
mw  = 20; % markersize

%%% Select example
%CAS = 'siso_passive_simple'; r = 1e-12; ds = 0;
CAS = 'siso_passive_aca'; r = 1e-12; ds = 0;
%CAS = 'siso_passive_gugercin'; r = 5; ds = 0;
%CAS = 'mimo_passive_aca'; r = 5; ds = 0
%CAS = 'msd'; r = 1e-12; ds=.1;
[G,S,nu,ny,ph]  = lf.examples(CAS);
[SZ_v,SZ,~,~]   = lf.spectral_zeros(S);

%%% Data
nip     = 100;
puls    = logspace(-1,2,nip);
la      = puls(1:2:end)*1i; la = sort([la conj(la)]);
mu      = puls(2:2:end)*1i; mu = sort([mu conj(mu)]);
% pts     = [1 2];
% la      = pts;
% mu      = -pts;
k       = length(la);
q       = length(mu);
R       = ones(nu,k);
L       = ones(q,ny);
for ii = 1:k; W(1:ny,1:nu,ii) = G(la(ii)); end
for ii = 1:q; V(1:ny,1:nu,ii) = G(mu(ii)); end

%%% Loewner 
opt                 = [];
opt.target          = r;
opt.real            = true;
[hloe,info_loe]     = lf.loewner_tng(la,mu,W,V,R,L,opt);
isPassive(info_loe.Hr)

%%% Loewner pH
D                   = S.D;
opt.Ds              = ds;
[hloep,info_loep]   = lf.loewner_passive(la,mu,W,V,R,L,D,opt);
[hloeph,info_loeph] = lf.passive2ph(info_loep.Hrn);

%%% Plot 
% >> Singular values
figure, hold on, grid on
plot(info_loe.sv,'-o','MarkerSize',mw,'LineWidth',lw)
plot(info_loe.sv_nu,'-o','MarkerSize',mw,'LineWidth',lw)
plot(info_loep.sv,'--x','MarkerSize',mw,'LineWidth',lw)
plot(info_loep.sv_nu,'--x','MarkerSize',mw,'LineWidth',lw)
set(gca,'YScale','log')
xlabel('$k$'), ylabel('Normalized singular value')
title(['Singular values $r=' num2str(r) '$'])
legend({'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)' ...
        'svd($[\bf{L},\bf{M}]$) pH','svd($\bf{L}$) pH'},'Location','best')
% >> Eigenvalues
eigS    = eig(S);
eigHt   = eig(info_loe.Hr);
eigHp   = eig(info_loep.Hr);
figure, hold on, grid on
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHt),imag(eigHt),'o','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHp),imag(eigHp),'s','MarkerSize',mw,'LineWidth',lw);
title(['Eigenvalues $r=' num2str(r) '$'])
xlabel('Real'), ylabel('Imag.')
legend({'Original' ['Loewner (\texttt{isPassive}:' num2str(isPassive(info_loe.Hr)) ')'] ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(info_loep.Hr)) ')']},'location','best')
% >> Bode magnitude
w = logspace(-2,3,1e4);
for i = 1:numel(w)
    G_(1:ny,1:nu,i)         = G(1i*w(i));
    hloe_(1:ny,1:nu,i)      = hloe(1i*w(i));
    hloep_(1:ny,1:nu,i)     = hloep(1i*w(i));
    hloeph_(1:ny,1:nu,i)    = hloeph(1i*w(i));
end
figure, hold on, kk = 0;
for out = 1:ny
    for in = 1:nu
        kk = kk + 1;
        subplot(ny,nu,kk), hold on, grid on
        plot(w,20*log10(abs(squeeze(G_(out,in,:)))),'-'), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hloe_(out,in,:)))),'--'), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hloep_(out,in,:)))),':'), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hloeph_(out,in,:)))),'-.'), set(gca,'XScale','log')
        if out==ny; xlabel('pulsation [rad/s]'); end
    end
end
sgtitle('Bode magnitude')
legend({'Original' ... 
       ['Loewner (\texttt{isPassive}:' num2str(isPassive(info_loe.Hr)) ')'] ...
       ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(info_loep.Hr)) ')'] ... 
       'pH Loewner'},'location','best')

info_loeph

%%% Time-domain simulation
% >> Projector
out     = 1;
in      = 1;
Hr      = info_loep.Hrn;
[Vproj,Wproj,Vproj_x0]  = lf.projectors(S,info_loep,'none');
% >> t, u, x0
f           = .01;
dt          = .01;
t           = 0:dt:20;
u           = (sin(2*pi*f*t.^3).*exp(-.1.*t))';
x0          = zeros(length(S.A),1);
% >> Original
[yy,tt,xx]  = lsim(S(out,in),u,t,x0);
x0r         = Vproj_x0*x0;
% >> Identified and lift
[yr,~,xr]   = lsim(Hr(out,in),u,t,x0r);
xrp         = (Vproj*xr.')';
% >> Compare 
figure, 
subplot(211), hold on
plot(tt,yy,'-','LineWidth',lw,'DisplayName','Original'), grid on
plot(tt,yr,'--','LineWidth',lw,'DisplayName','pH-ROM')
Lgnd = legend('show');
xlabel('$t$'), ylabel('Output'),% set(gca,'XScale','log')
subplot(212), hold on
h1=plot(tt,xx,'-','LineWidth',lw); grid on
h2=plot(tt,xr,'m:','LineWidth',lw);
h3=plot(tt,xrp,'k--','LineWidth',lw);
leg = legend([h1(1), h2(1), h3(1)], ... 
             'Original',...
             'pH-ROM', ...
             'pH-ROM (lifted)');
set(leg, 'interpreter','latex')
xlabel('$t$'), ylabel('Internal variables'),% set(gca,'XScale','log')

