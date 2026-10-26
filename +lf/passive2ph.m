% Passive to pH form algorithm
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
% 
% Syntax
% [hr,info] = lf.passive2ph(H)
%  
% Input arguments
%  - H : passive normalized state-space ODE (ss)
% 
% Output arguments
%  - hr   : approximation model (handle function)
%  - info : structure with informations about the Loewner world
%    * r    : rational order (integer)
%    * Hr   : state-space model (ss)
%    * Q    : energy matrix 
%    * J,N  : reflective matrices
%    * R,S  : resistive matrices
%    * G,P  : satisfy Wcal >= 0
%    * Wcal : [R P; P' S] >= 0
%    * Vcal : [-Jr -Gr; Gr' Nr] skew sym
%    * dx   : handle function of state equations
%    * y    : handle function of output equations
% 
% Description
% Loewner rules.
%

function [hloeph,info] = passive2ph(H)

%%% pH form in real struture
[A,B,C,D]   = ssdata(H);
Sblk        = [-A -B; C D];
n           = length(A);
JGN         = 1/2*(Sblk-Sblk.');
RPS         = 1/2*(Sblk+Sblk.');
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
info.r      = length(A);
info.Hr     = ss((Jr-Rr),(Gr-Pr),(Gr+Pr)',(Nr+Sr));
info.Q      = Qr;
info.J      = Jr;
info.R      = Rr;
info.G      = Gr;
info.P      = Pr;
info.N      = Nr;
info.S      = Sr;
info.Wcal   = Wcal;
info.Vcal   = Vcal;
%
info.dx     = @(t,x,u) ((info.J-info.R)*info.Q)*x + (info.G-info.P)*u;
info.y      = @(x,u) ((info.G+info.P).'*info.Q)*x + (info.N+info.S)*u;

