function [hloeph,info] = passive2ph(H)

%%% pH form in real struture
[A,B,C,D]   = ssdata(H);
Sblk        = [-A -B; C D];
n           = length(A);
JGN         = 1/2*(Sblk-Sblk');
RPS         = 1/2*(Sblk+Sblk');
Jr          = -JGN(1:n,1:n);
Gr          = -JGN(1:n,n+1:end);
Nr          = JGN(n+1:end,n+1:end);
Rr          = RPS(1:n,1:n);
Pr          = RPS(1:n,n+1:end);
Sr          = RPS(n+1:end,n+1:end);
Qr          = eye(n);
In          = eye(n);
hloeph      = @(s) ((Gr+Pr)')*((s*In-(Jr-Rr))\(Gr-Pr))+(Nr+Sr);
Wcal        = [ Rr  Pr; Pr' Sr];
Vcal        = [-Jr -Gr; Gr' Nr];
%
info.Q      = Qr;
info.J      = Jr;
info.R      = Rr;
info.G      = Gr;
info.P      = Pr;
info.N      = Nr;
info.S      = Sr;
info.Wcal   = Wcal;
info.Vcal   = Vcal;
