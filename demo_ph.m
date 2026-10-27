clearvars; close all; clc;
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_facto = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_facto,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_facto{index_interpreter(i)},'factory','default'),'latex'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BEGIN: TO ADAPT BY USER
% Add Loewner package
addpath('/Users/charles/Documents/GIT/lf')
%%% Choose case & order
%CAS     = 'heat';
%CAS     = 'siso_passive_gugercin'; robj = 5;
%CAS     = 'siso_passive_simple'; robj = 1e-12;
%CAS     = 'siso_passive_aca'; robj = []; %robj = 5
CAS     = 'msd'; robj = [];
%CAS     = 'mimo_passive_aca'; robj = [];
%%% END: TO ADAPT BY USER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% General variables 
lw      = 3;  % linewidth
mw      = 20; % markersize
colref  = orderedcolors('gem12');

%%% Select example
[G,S,nu,ny,ph]  = lf.examples(CAS);
[lf.isPositive(ph.Q) ... 
 lf.isPositive(ph.Wcal) ... 
 lf.isHermitian(ph.Wcal) ...
 lf.isSkewHermitian(ph.Vcal)]

%%% Data
w_  = 1i*logspace(-2,2,50);
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
opt.D       = 0*ones(ny,nu);
opt.real    = true;
[hloe,iloe] = lf.loewner_tng(la_,mu_,W,V,R,L,opt);

%%% Loewner pH
D               = S.D;
opt.Ds          = .1;
[hloep,iloep]   = lf.loewner_passive(la_,mu_,W,V,R,L,D,opt);
[hloeph,phr]    = lf.passive2ph(iloep.Hrn);
r               = length(phr.Q);
%
[lf.isPositive(phr.Q) ... 
 lf.isPositive(phr.Wcal) ... 
 lf.isHermitian(phr.Wcal) ...
 lf.isSkewHermitian(phr.Vcal)]

%%% Computed associated Petrov-Galerkin projectors
[Vp,Wp,Vx0] = lf.projectors(S,iloep,'none');
Wp          = Wp/((Wp'*Vp)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% >> Singular values
subplot(131), hold on, grid on
plot(iloe.sv,'-o','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($[\bf{L},\bf{M}]$)')
plot(iloe.sv_nu,'-x','MarkerSize',mw,'LineWidth',lw,'DisplayName','svd($\bf{L}$)')
plot(iloep.sv,'--^','MarkerSize',mw,'LineWidth',lw,'MarkerFaceColor',colref(7,:),'DisplayName','svd($[\bf{L},\bf{M}]$) passive')
plot(iloep.sv_nu,'--v','MarkerSize',mw,'LineWidth',lw,'MarkerFaceColor',colref(6,:),'DisplayName','svd($\bf{L}$) passive')
set(gca,'YScale','log')
xlabel('$k$','Interpreter','latex')
ylabel('Normalized singular value','Interpreter','latex')
title('Singular values','interpreter','latex','FontSize',26)
legend('show','interpreter','latex','location','best')
% >> Spectral zeros
subplot(132), hold on, grid on
[~,SZ]      = lf.spectral_zeros(S);
[~,SZloep]  = lf.spectral_zeros(iloep.Hr);
plot(real(SZ),imag(SZ),'.','MarkerSize',mw,'LineWidth',lw,'DisplayName','FOM');
plot(real(SZloep),imag(SZloep),'o','MarkerSize',mw,'LineWidth',lw,'DisplayName','ROM (passive)');
axis tight
xlabel('Real part','Interpreter','latex')
ylabel('Imag. part','Interpreter','latex')
legend('show','Interpreter','latex','location','best')
title('Spectral zeros','interpreter','latex','FontSize',26)
% >> Eigenvalues
subplot(133), hold on, grid on
eigS        = eig(S);
eigHloep    = eig(iloep.Hr);
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw,'DisplayName','FOM');
plot(real(eigHloep),imag(eigHloep),'o','MarkerSize',mw,'LineWidth',lw,'DisplayName','ROM (passive)');
xlabel('Real part','Interpreter','latex')
ylabel('Imag. part','Interpreter','latex')
legend('show','Interpreter','latex','location','best')
title('Eigenvalues','interpreter','latex','FontSize',26)
% >> Surface in the SISO case only
if nu == 1
    lim_xx  = max(abs(union(la_,mu_)));
    xSpace  = linspace(-lim_xx*1.1,lim_xx*1.1,2e2);
    ySpace  = linspace(-lim_xx*1.1,lim_xx*1.1,2e2);
    for i = 1:numel(xSpace)
        for j = 1:numel(ySpace)
            tab(i,j)  = G(xSpace(i)+1i*ySpace(j));
            tabr(i,j) = hloeph(xSpace(i)+1i*ySpace(j));
        end
    end
    [X,Y] = meshgrid(xSpace,ySpace);
    figure, hold on, grid on,
    surf(X,Y,log10(abs(tab.')),'FaceColor','none');%,'DisplayName','$\mathbf H$')
    surfc(X,Y,log10(abs(tabr.')),'EdgeColor','none');%,'DisplayName','$\mathbf{H}$')
    plot3(real(la_),imag(la_),log10(abs(W(:))),'s','MarkerSize',mw,'MarkerFaceColor','k','LineWidth',lw)
    plot3(real(mu_),imag(mu_),log10(abs(V(:))),'o','MarkerSize',mw,'MarkerFaceColor','k','LineWidth',lw)
    h = gca();
    plot3(real(eigS),imag(eigS),ones(1,length(eigS))*h.ZLim(2),'d','MarkerSize',mw,'LineWidth',lw)%,'MarkerFaceColor','k')
    plot3(real(eigHloep),imag(eigHloep),ones(1,length(eigHloep))*h.ZLim(2),'+','MarkerSize',mw,'LineWidth',lw)%,'MarkerFaceColor','k')
    legend({'$\mathbf H$','$\mathbf{\hat{H}}$ (pH)','$\mathbf{\hat{H}}$ (pH) contour','$\{\lambda_j,\mathbf w_j\}$','$\{\mu_i,\mathbf v_i\}$','$\mathbf{H}$ eigen.','$\mathbf{\hat{H}}$ eigen.'})
    view(20,35), alpha(.5)
    xlabel('Real part $x$'); ylabel('Imag. part $y$'); zlabel('Magnitude dB, $\log|\mathbf H(x+\imath y)|$'); 
    title('Function evaluation over the complex domain','interpreter','latex','FontSize',26)
end
% >> Bodemag
figure, hold on, kk = 0;
w = logspace(-2,3,300);
for i = 1:numel(w)
    G_(1:ny,1:nu,i)         = G(1i*w(i));
    hloep_(1:ny,1:nu,i)     = hloep(1i*w(i));
    hloeph_(1:ny,1:nu,i)    = hloeph(1i*w(i));
end
for out = 1:ny
    for in = 1:nu
        kk = kk + 1;
        subplot(ny,nu,kk), hold on, grid on
        plot(w,20*log10(abs(squeeze(G_(out,in,:)))),'-','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hloep_(out,in,:)))),'--','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        plot(w,20*log10(abs(squeeze(hloeph_(out,in,:)))),':','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        if out==ny; xlabel('Pulsation [rad/s]'); end
        if in==1; ylabel('Gain [dB]'); end
    end
end
sgtitle('Bode magnitude','interpreter','latex','FontSize',26)
legend({'pH-FOM' 'ROM (passive)' 'pH-ROM'},'Interpreter','latex','location','best')
% >> Nyquist
figure, hold on, kk = 0;
for out = 1:ny
    for in = 1:nu
        kk = kk + 1;
        subplot(ny,nu,kk), hold on, grid on
        plot(real(squeeze(G_(out,in,:)))     ,imag(squeeze(G_(out,in,:))),'-','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        plot(real(squeeze(hloep_(out,in,:))) ,imag(squeeze(hloep_(out,in,:))),'--','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        plot(real(squeeze(hloeph_(out,in,:))),imag(squeeze(hloeph_(out,in,:))),':','MarkerSize',mw,'LineWidth',lw), set(gca,'XScale','log')
        if out==ny; xlabel('real'); end
        if in==1; ylabel('imag'); end
    end
end
sgtitle('Nyquist','interpreter','latex','FontSize',26)
legend({'pH-FOM' 'ROM (passive)' 'pH-ROM'},'Interpreter','latex','location','best')
% >> Time-domain simulation
f           = .1;
dt          = .01;
t           = 0:dt:15;
u           = @(t) (sin(2*pi*f*t.^3).*exp(-.5*t)).*ones(nu,1);
uu          = u(t);
x0          = zeros(length(S.A),1); %x0(1)=1e-2;
% Original
[tt,xx]     = ode45(@(t,x) ph.dx(t,x,u(t)), t, x0); xx = xx'; tt = tt';
yy          = ph.y(xx,uu);
% Identified and lifted
x0r         = Vx0*x0;
[tt,xr]     = ode45(@(t,x) phr.dx(t,x,u(t)), tt, x0r); xr = xr'; tt = tt';
yr          = phr.y(xr,uu);
xrp         = Vp*xr;
% Compare output and internal variables
figure, 
subplot(311), hold on
plot(tt,uu,'-','LineWidth',lw,'DisplayName','$u$'), grid on
legend('show'); xlabel('Time [s]'), ylabel('Input'),
subplot(312), hold on
plot(tt,yy,'-','Color',colref(1,:),'LineWidth',lw,'DisplayName','pH-FOM'), grid on
plot(tt,yr,'--','Color',colref(2,:),'LineWidth',lw,'DisplayName','pH-ROM')
legend('show'); xlabel('Time [s]'), ylabel('Output'),
subplot(313), hold on
h1=plot(tt,xx,'-','LineWidth',lw,'Color',colref(1,:)); grid on
h2=plot(tt,xr,'--','Color',colref(2,:),'LineWidth',lw);
h3=plot(tt,xrp,':','Color',colref(3,:),'LineWidth',lw);
leg = legend([h1(1), h2(1), h3(1)], ... 
             ['pH-FOM, $\mathbf x\in R^{' num2str(size(xx,1)) '}$'],...
             ['pH-ROM, $\mathbf{\hat x}\in R^{' num2str(size(xr,1)) '}$'],...
             ['pH-ROM (lifted), $\mathbf{\hat x}_{\uparrow}=V\mathbf{\hat x}\in R^{' num2str(size(xrp,1)) '}$']);
set(leg, 'interpreter','latex')
xlabel('Time [s]'), ylabel('Internal variables'),% set(gca,'XScale','log')
sgtitle('Time-domain response','interpreter','latex','FontSize',26)
%
n = size(xx,1);
if n < 10
    figure
    for ii = 1:n
        subplot(ceil(n/2),2,ii), hold on
        plot(tt,abs(xx(ii,:)),'-','LineWidth',lw,'Color',colref(1,:),'DisplayName',['pH-FOM, $\mathbf x\in R^{' num2str(size(xx,1)) '}$']); grid on
        plot(tt,abs(xrp(ii,:)),'--','LineWidth',lw,'Color',colref(2,:),'DisplayName',['pH-ROM (lifted), $\mathbf{\hat x}_{\uparrow}=V\mathbf{\hat x}\in R^{' num2str(size(xrp,1)) '}$ and $\mathbf{\hat x}\in R^{' num2str(size(xr,1)) '}$']);
        plot(tt,abs(xx(ii,:)-xrp(ii,:)),'k-','LineWidth',lw,'DisplayName','error'); grid on, set(gca,'YScale','log')
        xlabel('Time [s]')
        h=gca; set(gca,'YLim',[1e-9 max(h.YLim)])
    end
    legend('show','location','best')
    sgtitle('Time-domain internal variables responses','interpreter','latex','FontSize',26)
end
% Energy
for ii = 1:numel(tt)
    %
    Efom(ii)    = uu(:,ii)'*yy(:,ii) - uu(:,ii)'*ph.S*uu(:,ii);
    Erom(ii)    = uu(:,ii)'*yr(:,ii) - uu(:,ii)'*phr.S*uu(:,ii);
    % 
    Hfom(ii)    = .5*xx(:,ii)'*ph.Q*xx(:,ii);
    Hrom(ii)    = .5*xr(:,ii)'*phr.Q*xr(:,ii);
    Hromp(ii)   = .5*xrp(:,ii)'*(Vp*phr.Q*Wp')*xrp(:,ii);
    %
    Dfom(ii)    = -xx(:,ii)'*ph.R*xx(:,ii);% + 2*xx(:,ii)'*ph.Q*ph.P*uu(:,ii);
    Drom(ii)    = -xr(:,ii)'*phr.R*xr(:,ii);% + 2*xr(:,ii)'*phr.Q*phr.P*uu(:,ii);
    Dromp(ii)   = -xrp(:,ii)'*(Vp*phr.R*Wp')*xrp(:,ii);
    %
    Hfom_dt(ii) = -xx(:,ii)'*ph.Q'*ph.R*ph.Q*xx(:,ii)                 + uu(:,ii)'*yy(:,ii) - uu(:,ii)'*ph.S *uu(:,ii);% - 2*xx(:,ii)'*ph.Q*ph.P*uu(:,ii);
    Hrom_dt(ii) = -xr(:,ii)'*phr.Q'*phr.R*phr.Q*xr(:,ii)              + uu(:,ii)'*yr(:,ii) - uu(:,ii)'*phr.S*uu(:,ii);% - 2*xr(:,ii)'*phr.Q*phr.P*uu(:,ii);
    Hromp_dt(ii)= -xrp(:,ii)'*(Vp*(phr.Q'*phr.R*phr.Q)*Wp')*xrp(:,ii) + uu(:,ii)'*yr(:,ii) - uu(:,ii)'*phr.S*uu(:,ii);% - 2*xr(:,ii)'*phr.Q*phr.P*uu(:,ii);
end
%
figure, 
subplot(221), hold on
plot(tt,Hfom ,'-','LineWidth',lw,'DisplayName','$\mathcal H=\frac{1}{2}\mathbf x^\top Q \mathbf x$'), grid on
plot(tt,Hrom ,'--','LineWidth',lw,'DisplayName','$\hat \mathcal H=\frac{1}{2}\mathbf{\hat{x}}^\top \hat Q \mathbf{\hat{x}}$'),
plot(tt,Hromp,':','LineWidth',lw,'DisplayName','$\hat \mathcal H_{\uparrow}=\frac{1}{2}\mathbf{\hat{x}}_{\uparrow}^\top V \hat Q W^\top\mathbf{\hat{x}}_{\uparrow}$'),
legend('show','Location','Best'); xlabel('Time [s]'), ylabel('Hamiltonian'),
subplot(222), hold on
plot(tt,Efom,'-','LineWidth',lw,'DisplayName','$\mathcal E=\mathbf{u}^\top\mathbf{y}-\mathbf{u}^\top S \mathbf{u}$'), grid on
plot(tt,Erom,'--','LineWidth',lw,'DisplayName','$\hat \mathcal E=\mathbf{\hat u}^\top\mathbf{\hat y}-\mathbf{u}^\top \hat S \mathbf{u}$')
legend('show','Location','Best'); xlabel('Time [s]'), ylabel('Energy supplier'),
subplot(223), hold on
plot(tt,Dfom ,'-','LineWidth',lw,'DisplayName','$\mathcal D=-\mathbf x^\top R \mathbf x$'), grid on
plot(tt,Drom ,'--','LineWidth',lw,'DisplayName','$\hat \mathcal D=-\mathbf{\hat{x}}^\top \hat R \mathbf{\hat{x}}$')
plot(tt,Dromp,':','LineWidth',lw,'DisplayName','$\hat \mathcal D=-\mathbf{\hat{x}}_{\uparrow}^\top V \hat R W^\top \mathbf{\hat{x}_{\uparrow}}$')
legend('show','Location','Best'); xlabel('Time [s]'), ylabel('Damping'),
subplot(224), hold on
plot(tt,Efom+Dfom -Hfom_dt ,'-','LineWidth',lw,'DisplayName','$ \mathcal E+ \mathcal D - \frac{d}{dt}\mathcal H$'), grid on
plot(tt,Erom+Drom -Hrom_dt ,'--','LineWidth',lw,'DisplayName','$ \hat \mathcal E+ \hat \mathcal D -\frac{d}{dt}\hat\mathcal H$')
plot(tt,Erom+Dromp-Hromp_dt,':','LineWidth',lw,'DisplayName','$ \hat \mathcal E_{\uparrow}+ \hat \mathcal D_{\uparrow} -\frac{d}{dt}\hat\mathcal H_{\uparrow}$')
legend('show','Location','Best'); xlabel('Time [s]'), ylabel('Balance'),
sgtitle(['Energy balance ($n=' num2str(n) '$ and $r=' num2str(r) '$)'],'interpreter','latex','FontSize',26)
