function [G,S,nu,ny,eigS] = examples(CAS)

switch CAS
    case 'siso_vsimple' % Very simple
        ny  = 1;
        nu  = 1;
        S   = ss(tf([1 1],[1 5]));
        S.E = eye(length(S.A));
    case 'siso_simple' % simple
        ny  = 1;
        nu  = 1;
        S   = ss(tf(1,[1 0.01 1]));
        S.E = eye(length(S.A));
    case 'mimo_rand'
        ny  = 4;
        nu  = 3;
        S   = rss(20,ny,nu);
        S.E = eye(length(S.A));
        S.D = 2*ones(ny,nu);
    case 'mimo_large' % Full I/O 
        n   = 7;
        nu  = n;
        ny  = n;
        S   = rss(n,ny,nu);
        S.E = eye(n);
        S.B = eye(n);
        S.C = eye(n);
        S.D = ones(n,n);
    %%% PASSIVE
    case 'siso_passive_simple' % Very simple
        ny  = 1;
        nu  = 1;
        S   = ss(tf([2 4],[1 1]));
        S.E = eye(length(S.A));
    case 'siso_passive_aca' % Passive RLC by ACA
        A   = [-20 -10 0 0 0; 10 0 -10 0 0; 0 10 0 -10 0; 0 0 10 0 -10; 0 0 0 10 -2];
        B   = [20 0 0 0 0]';
        nu  = size(B,2);
        ny  = nu;
        C   = B';
        D   = 2*eye(nu);
        E   = eye(size(A));
        S   = dss(A,B,C,D,E);
    case 'mimo_passive_aca' % Passive RLC by ACA
        A   = [-20 -10 0 0 0; 10 0 -10 0 0; 0 10 0 -10 0; 0 0 10 0 -10; 0 0 0 10 -2];
        B   = [20 0 0 0 0; 0 0 20 0 0; 0 0 0 0 1]';
        nu  = size(B,2);
        ny  = nu;
        C   = B';
        D   = 2*eye(nu);
        E   = eye(size(A));
        S   = dss(A,B,C,D,E);
    case 'siso_passive_gugercin' % Passive RLC by Serkan
        load('+lf/rlc_serkan200.mat')
        nu  = size(B,2);
        ny  = nu;
        E   = eye(size(A));
        S   = dss(A,B,C,D,E);
end
%
eigS        = eig(S);
[A,B,C,D,E] = dssdata(S);
G           = @(s) C*((s*E-A)\B)+D;