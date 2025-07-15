function [uc,tc,info] = multisine(Ns,Ts,FBND,REV,SHOW)

% Nyquist check
Fs  = 1/Ts;
if max(FBND) >= Fs/2
    error('Nyquist frequency issue: check "max(fBnd) < 1/(2Ts)".')
end
%
tc  = (0:Ns-1)*Ts;
f0  = Fs/Ns;
% Select bounds
fId = 1 + round(FBND./f0);
%F   = length(fId(1):fId(2));
F   = fId(2) - fId(1) + 1;
amp = zeros(1,fId(2));
amp(end-F+1:end) = 1;
% 
u   = zeros(F,Ns);
for k = 1:length(amp)
    %phi_k   = rand(1)*2*pi; % Random multisine
    phi_k   = -k*(k-1)*pi/F; % Schroeder multisine
    u(k,:)  = amp(k)*1/F*cos(2*pi*f0*k*tc + phi_k);
end
uc   = sum(u,1);
uc   = uc/max(uc(:));
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
%
info.uk     = u;
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
    legend({'Continuous-time'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    %
    subplot(212); hold on; grid on, axis tight
    plot(fc,abs(FTuc),'LineWidth',3),
    hh = gca;
    plot([1 1]*Fs/2,[hh.YLim(1) hh.YLim(2)],'k:','LineWidth',3), 
    set(gca,'TickLabelInterpreter','latex','FontSize',FONT_SZ2)
    xlabel('$f$ [Hz]','Interpreter','latex','FontSize',FONT_SZ), 
    ylabel('$\mathbf{U}(f)$','Interpreter','latex','FontSize',FONT_SZ)
    legend({'Continuous-time (FFT)','Nyquist frequency'},'Location','East','Interpreter','latex','FontSize',FONT_SZ)
    %
    sgtitle(['Multisine signal $\{N_s,T_s,T_f\}=\{' num2str(Ns) ',' num2str(Ts)  ',' num2str(tc(end)) '\}$'],'Interpreter','latex','Fontsize',20)
end
