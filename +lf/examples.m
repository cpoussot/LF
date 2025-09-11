function [G,S,nu,ny,eigS] = examples(CAS)

switch CAS
    case 'mimo_rand'
        ny  = 4;
        nu  = 3;
        S   = rss(20,ny,nu);
        S.E = eye(length(S.A));
        S.D = 2*ones(ny,nu);
    case 'siso_passive' % Passive RLC
        A   = [-20 -10 0 0 0; 10 0 -10 0 0; 0 10 0 -10 0; 0 0 10 0 -10; 0 0 0 10 -2];
        B   = [20 0 0 0 0]';
        nu  = size(B,2);
        ny  = nu;
        C   = B';
        D   = 1*eye(nu);
        E   = eye(size(A));
        S   = dss(A,B,C,D,E);
    case 'siso_simple' % Very simple
        ny  = 1;
        nu  = 1;
        S   = ss(tf(1,[1 0.01 1]));
        S.E = eye(length(S.A));
        nip = 4;
    case 'mimo_large' % Full I/O 
        n   = 7;
        nu  = n;
        ny  = n;
        S   = rss(n,ny,nu);
        S.E = eye(n);
        S.B = eye(n);
        S.C = eye(n);
        S.D = ones(n,n);
end
%
eigS        = eig(S);
[A,B,C,D,E] = dssdata(S);
G           = @(s) C*((s*E-A)\B)+D;