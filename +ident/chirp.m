% <strong>INSA DYNAMICAL SYSTEMS SUITE</strong>
% Main interface for chirp signal generation
% Author: C. Poussot-Vassal [MOR Digital Systems / Onera]
% Date  : April 2021 (creation)
%
% Description
% Computes a chirp signal sweeping from f0 to f1.
%
% Syntax
%  u_chirp = insa.chirp(t,f0,t1,f1)
%  u_chirp = insa.chirp(t,f0,t1,f1,type)
% 
% Input arguments
%  - t  : time vector with constant sampling Ts (real, N x 1)
%  - f0 : starting frequency in Hz (real)
%  - t1 : final frequency time (real <= t(end))
%  - f1 : final frequency in Hz (real)
% 
% Optional arguments
%  - type : variation of the chirp linear from f0 -> f1 being 'linear'
%           (default), 'quadratic' or 'logarithmic' 
% 
% Output arguments
%  - u : chirp signal in [-1,1]
% 

function [uc,tc,info] = chirp(Ns,Ts,FBND,REV,TYPE,SHOW)

% Nyquist check
Fs  = 1/Ts;
if max(FBND) >= Fs/2
    error('Nyquist frequency issue: check "max(fBnd) < 1/(2Ts)".')
end
%
tc  = (0:Ns-1)*Ts;
f0  = min(FBND);
f1  = max(FBND);
t1  = tc(end);
% Set p=1 for linear, 2 for quadratic, 3 for logarithmic
switch TYPE
    case 'linear'
        p = 1;
    case 'quadratic'
        p = 2;
    case 'logarithmic'
        p = 3;
    otherwise
        p = 1;
end
beta    = (f1-f0).*(t1.^(-p));
uc      = sin(2*pi * ( beta./(1+p).*(tc.^(1+p)) + f0.*tc));
if REV
    uc = fliplr(uc);
end
% FFT
L       = numel(uc);
FTuc    = fft(uc)/L;
fc      = linspace(0,1,L)*Fs;
% Info
info.Ts     = Ts;
info.Fs     = Fs;
info.fc     = fc;
info.FTuc   = FTuc;
% Plot
if SHOW
    FONT_SZ     = 16;
    FONT_SZ2    = 14;
    %
    figure, 
    subplot(211); hold on, grid on, axis tight
    plot(tc,uc,'-','LineWidth',3),
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$t$ [s]','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{u}(t)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Time-domain'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    %
    subplot(212); hold on; grid on, axis tight
    plot(fc,abs(FTuc),'LineWidth',3),
    hh = gca;
    plot([1 1]*Fs/2,[hh.YLim(1) hh.YLim(2)],'k:','LineWidth',3), 
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$f$ [Hz]','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{U}(f)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'FFT','Nyquist frequency'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    %
    sgtitle(['Frequency chirp signal $\{N_s,T_s,T_f\}=\{' num2str(Ns) ',' num2str(Ts)  ',' num2str(tc(end)) '\}$'],'Interpreter','latex','Fontsize',20)
end
