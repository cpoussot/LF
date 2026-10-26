clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Loewner package
addpath('/Users/charles/Documents/GIT/lf')
%%% Select example
CAS = 'siso_simple';
CAS = 'mimo_rand';
%CAS = 'mimo_large';
%%% Select order
robj = 1e-12;
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% General variables 
lw  = 3;  % linewidth
mw  = 20; % markersize

%%% Select example
[G,S,nu,ny] = lf.examples(CAS);
size(S)

%%% Data
nip = 100;
la_ = (logspace(-2,2,nip))*1i;    la_ = sort([la_ conj(la_)]);
mu_ = (logspace(-2,2,nip)+.1)*1i; mu_ = sort([mu_ conj(mu_)]);
k   = length(la_);
q   = length(mu_);
%R   = ones(nu,k);
%L   = ones(q,ny);
for ii = 1:k; W(1:ny,1:nu,ii) = G(la_(ii)); end
for ii = 1:q; V(1:ny,1:nu,ii) = G(mu_(ii)); end
for ii = 1:2:k; R(1:nu,ii) = rand(nu,1)+1i*rand(nu,1); R(1:nu,ii+1) = conj(R(1:nu,ii)); end
for ii = 1:2:q; L(ii,1:ny) = rand(1,ny)+1i*rand(1,ny); L(ii+1,1:ny) = conj(L(ii,1:ny)); end

%%% Loewner 
opt         = [];
opt.target  = robj;
opt.D       = 0*ones(ny,nu);
opt.real    = true;
% >> Tangential 
[htng,itng] = lf.loewner_tng(la_,mu_,W,V,R,L,opt);
if ~itng.isCC
    test1 = itng.MU*itng.LL - itng.LL*itng.LA;
    test2 = itng.V*itng.R - itng.L*itng.W;
    norm(test1-test2)
    test1 = itng.MU*itng.SS - itng.SS*itng.LA;
    test2 = itng.MU*itng.V*itng.R - itng.L*itng.W*itng.LA;
    norm(test1-test2)
end
% >> Block
[hblk,iblk] = lf.loewner_blk(la_,mu_,W,V,opt);

%%% Some informations
itng
iblk
size(itng.Hr)
size(iblk.Hr) 
[isreal(itng.Hr) isreal(iblk.Hr)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >> Singular values
figure, hold on, grid on
plot(itng.sv,'-o','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($[\bf{L},\bf{M}]$)')
plot(itng.sv_nu,'-x','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($\bf{L}$)')
plot(iblk.sv,'--o','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($[\bf{L},\bf{M}]$)')
plot(iblk.sv_nu,'--x','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($\bf{L}$)')
set(gca,'YScale','log')
xlabel('$k$','Interpreter','latex')
ylabel('Normalized singular value','Interpreter','latex')
title('Normalized singular values')
legend('show','location','best','interpreter','latex')
% >> Eigenvalues
eigS  = eig(S);
eigHt = eig(itng.Hr);
eigHb = eig(iblk.Hr);
figure, hold on, grid on
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw,'DisplayName','Original');
plot(real(eigHt),imag(eigHt),'o','MarkerSize',mw,'LineWidth',lw,'DisplayName','Loewner tangent');
plot(real(eigHb),imag(eigHb),'s','MarkerSize',mw,'LineWidth',lw,'DisplayName','Loewner block' );
xlabel('Real','Interpreter','latex')
ylabel('Imag.','Interpreter','latex')
title('Eigenvalues')
legend('show','location','best','interpreter','latex')
% >> Bodemag
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
sgtitle('Bode magnitude','FontSize',24)
legend({'Original' 'Loewner tangent' 'Loewner block'},'Interpreter','latex','location','best')
