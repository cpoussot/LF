clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_facto = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_facto,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_facto{index_interpreter(i)},'factory','default'),'latex'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Loewner package
addpath('D:\cpoussot\Documents\ONERA\GIT\LF')
%%% Choose case & order
%CAS = 'siso_simple'; robj = []
%CAS = 'siso_rand'; robj = []
%CAS = 'mimo_rand'; robj = 10;%[]
CAS = 'mimo_large'; robj = 10;
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% General variables 
lw      = 3;  % linewidth
mw      = 20; % markersize
colref  = orderedcolors('gem12');

%%% Select example
[G,S,nu,ny,mdl] = lf.examples(CAS);

%%% Data
w_  = 1i*logspace(-2,2,100);
la_ = w_(2:2:end); la_ = cplxpair([la_ conj(la_)]);
mu_ = w_(1:2:end); mu_ = cplxpair([mu_ conj(mu_)]);
k   = length(la_);
q   = length(mu_);
for ii = 1:k;   W(1:ny,1:nu,ii) = G(la_(ii)); end
for ii = 1:q;   V(1:ny,1:nu,ii) = G(mu_(ii)); end
for ii = 1:2:k; R(1:nu,ii) = rand(nu,1)+1i*rand(nu,1); R(1:nu,ii+1) = conj(R(1:nu,ii)); end
for ii = 1:2:q; L(ii,1:ny) = rand(1,ny)+1i*rand(1,ny); L(ii+1,1:ny) = conj(L(ii,1:ny)); end

%%% Loewner tangential
opt         = [];
opt.target  = robj;
opt.real    = true;
[hloe,iloe] = lf.loewner_tng(la_,mu_,W,V,R,L,opt);
r           = length(iloe.r);
mdlr.mass   = iloe.Er;
mdlr.dx     = @(t,x,u) iloe.Ar*x + iloe.Br*u;
mdlr.y      = @(x,u)   iloe.Cr*x + iloe.Dr*u;

%%% Computed associated Petrov-Galerkin projectors
[Vp,Wp,Vx0] = lf.projectors(S,iloe,'svd');
Wp          = Wp/((Wp'*Vp)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% >> Singular values
subplot(121), hold on, grid on
plot(iloe.sv,'-o','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($[\bf{L},\bf{M}]$)')
plot(iloe.sv_nu,'-x','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($\bf{L}$)')
set(gca,'YScale','log')
xlabel('$k$','Interpreter','latex')
ylabel('Normalized singular value','Interpreter','latex')
title('Singular values','interpreter','latex','FontSize',26)
legend('show','interpreter','latex','location','best')
% >> Eigenvalues
subplot(122), hold on, grid on
eigS    = eig(S);
eigHr   = eig(iloe.Hr);
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw,'DisplayName','FOM');
plot(real(eigHr),imag(eigHr),'o','MarkerSize',mw,'LineWidth',lw,'DisplayName','ROM');
xlabel('Real part','Interpreter','latex')
ylabel('Imag. part','Interpreter','latex')
legend('show','Interpreter','latex','location','best')
title('Eigenvalues','interpreter','latex','FontSize',26)
% >> Bodemag
figure, hold on, kk = 0;
w = logspace(-2,2,500);
for i = 1:numel(w)
    G_(1:ny,1:nu,i)     = G(1i*w(i));
    hloe_(1:ny,1:nu,i)  = hloe(1i*w(i));
end
for out = 1:ny
    for in = 1:nu
        kk = kk + 1;
        subplot(ny,nu,kk), hold on, grid on, axis tight
        plot(w,20*log10(abs(squeeze(G_(out,in,:)))),'-','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hloe_(out,in,:)))),'--','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        if out==ny; xlabel('Pulsation [rad/s]'); end
        if in==1; ylabel('Gain [dB]'); end
    end
end
sgtitle('Bode magnitude','interpreter','latex','FontSize',26)
legend({'FOM' 'ROM'},'Interpreter','latex','location','best')

%%
% >> Time-domain simulation
f           = .1;
dt          = .01;
t           = 0:dt:15;
u           = @(t) (sin(2*pi*f*t.^3).*exp(-.5*t)).*ones(nu,1);
uu          = u(t);
x0          = zeros(length(S.A),1); %x0(1)=1e-2;
% Original
[tt,xx]     = ode45(@(t,x) mdl.dx(t,x,u(t)), t, x0); xx = xx'; tt = tt';
yy          = mdl.y(xx,uu);
% Identified and lifted
x0r         = Vx0*x0;
options     = odeset('Mass',mdlr.mass);
[tt,xr]     = ode45(@(t,x) mdlr.dx(t,x,u(t)), tt, x0r, options); xr = xr'; tt = tt';
yr          = mdlr.y(xr,uu);
xrp         = Vp*xr;
% Compare output and internal variables
figure, 
subplot(311), hold on, grid on
plot(tt,uu,'-','LineWidth',lw,'DisplayName','$u$'), grid on
legend('show'); xlabel('Time [s]'), ylabel('Input'),
subplot(312), hold on, grid on
plot(tt,yy,'-','Color',colref(1,:),'LineWidth',lw,'DisplayName','FOM'), grid on
plot(tt,yr,'--','Color',colref(2,:),'LineWidth',lw,'DisplayName','ROM')
legend('show'); xlabel('Time [s]'), ylabel('Output'),
subplot(313), hold on, grid on
h1=plot(tt,xx,'-','LineWidth',lw,'Color',colref(1,:)); grid on
h2=plot(tt,xr,'--','Color',[1 1 1]*.9,'LineWidth',lw);
h3=plot(tt,xrp,':','Color',colref(3,:),'LineWidth',lw);
leg = legend([h1(1), h2(1), h3(1)], ... 
             ['FOM, $\mathbf x\in R^{' num2str(size(xx,1)) '}$'],...
             ['ROM, $\mathbf{\hat x}\in R^{' num2str(size(xr,1)) '}$'],...
             ['ROM (lifted), $\mathbf{\hat x}_{\uparrow}=V\mathbf{\hat x}\in R^{' num2str(size(xrp,1)) '}$']);
set(leg, 'interpreter','latex')
set(gca,'YLim',[min(xx(:)) max(xx(:))])
xlabel('Time [s]'), ylabel('Internal variables'),% set(gca,'XScale','log')
sgtitle('Time-domain response','interpreter','latex','FontSize',26)
%
n = size(xx,1);
if n < 10
    figure
    for ii = 1:n
        subplot(ceil(n/2),2,ii), hold on
        plot(tt,abs(xx(ii,:)),'-','LineWidth',lw,'Color',colref(1,:),'DisplayName',['FOM, $\mathbf x\in R^{' num2str(size(xx,1)) '}$']); grid on
        plot(tt,abs(xrp(ii,:)),'--','LineWidth',lw,'Color',colref(2,:),'DisplayName',['ROM (lifted), $\mathbf{\hat x}_{\uparrow}=V\mathbf{\hat x}\in R^{' num2str(size(xrp,1)) '}$ and $\mathbf{\hat x}\in R^{' num2str(size(xr,1)) '}$']);
        plot(tt,abs(xx(ii,:)-xrp(ii,:)),'k-','LineWidth',lw,'DisplayName','error'); grid on, set(gca,'YScale','log')
        xlabel('Time [s]')
        h=gca; set(gca,'YLim',[1e-9 max(h.YLim)])
    end
    legend('show','location','best')
    sgtitle('Time-domain internal variables responses','interpreter','latex','FontSize',26)
end
