clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end

%%% General variables 
lw  = 3;  % linewidth
mw  = 20; % markersize

%%% Select example
CAS = 'siso_simple';
%CAS = 'siso_passive_simple';
%CAS = 'mimo_rand';
%CAS = 'mimo_large';
[G,S,nu,ny] = lf.examples(CAS);

%%% Data
nip = 100;
la_ = (logspace(-2,2,nip))*1i;    la_ = sort([la_ conj(la_)]);
mu_ = (logspace(-2,2,nip)+.1)*1i; mu_ = sort([mu_ conj(mu_)]);
%la_ = (logspace(-2,2,nip));
%mu_ = [(logspace(-2,2,nip-2)+.1) -1i 1i];
k   = length(la_);
q   = length(mu_);
R   = ones(nu,k);
L   = ones(q,ny);
for ii = 1:k; W(1:ny,1:nu,ii) = G(la_(ii)); end
for ii = 1:q; V(1:ny,1:nu,ii) = G(mu_(ii)); end

%%% Loewner 
opt         = [];
opt.target  = 1e-12;
opt.D       = 0*ones(ny,nu);
opt.real    = true;
% Tangential 
[htng,itng] = lf.loewner_tng(la_,mu_,W,V,R,L,opt);
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

%%% Some informations
itng
iblk
size(itng.Hr)
size(iblk.Hr) 
[isreal(itng.Hr) isreal(iblk.Hr)]

%%% Plot
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

eigS  = eig(S);
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
for i = 1:numel(w)
    G_(1:ny,1:nu,i)     = G(1i*w(i));
    htng_(1:ny,1:nu,i)  = htng(1i*w(i));
    hblk_(1:ny,1:nu,i)  = hblk(1i*w(i));
end

figure, hold on, kk = 0;
for out = 1:ny
    for in = 1:nu
        kk = kk + 1;
        subplot(ny,nu,kk), hold on, grid on
        plot(w,20*log10(abs(squeeze(G_(out,in,:)))),'-'), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(htng_(out,in,:)))),'--'), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hblk_(out,in,:)))),':'), set(gca,'XScale','log')
        if out==ny; xlabel('pulsation [rad/s]'); end
    end
end
sgtitle('Bode magnitude')
legend({'Original' 'Loewner tangent' 'Loewner block'},'Interpreter','latex','location','best')
