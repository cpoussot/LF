function [G,S,nu,ny,eigS,ph] = examples(CAS)

ph = [];
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
    case 'siso_rand'
        rng(1234)
        ny  = 1;
        nu  = 1;
        S   = rss(100,ny,nu);
        S.E = eye(length(S.A));
        S.D = zeros(ny,nu);
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
    case 'msd'
        m = 10; b = 3; k = 1;
        m = 1; b = 1; k = 1;
        %
        n   = 2;
        ny  = 1;
        nu  = ny;
        % SS X=[x dx]
        E   = eye(n);
        A   = [0 1; -k/m -b/m];
        B   = [0; 1/m];
        C   = [0 1];
        D   = 0;
        S   = dss(A,B,C,D,E);
        % 
        In      = eye(n);
        ph.J    = [0 1; -1 0]/m;
        ph.R    = [0 0; 0 b]/m^2;
        ph.Q    = [k 0; 0 m];
        ph.G    = [0; 1/m];
        ph.P    = [0; 0];
        ph.N    = 0;
        ph.S    = 0;
        ph.h    = @(s) ((ph.G+ph.P).'*ph.Q)*((s*In-(ph.J-ph.R)*ph.Q)\(ph.G-ph.P))+(ph.N+ph.S);
end
%
eigS        = eig(S);
[A,B,C,D,E] = dssdata(S);
G           = @(s) C*((s*E-A)\B)+D;