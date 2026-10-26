function [A, B, C, D] = circuit_sparse(N, R , Rbar, Cap, L)
% circuit_sparse generates a SISO RLC ladder model of prescribed dimensions
%
% The model follows the descriptioin in DOI:10.1080/00207170410001713448
% and the input parameters R, Rbar, Cap and L default to the values given
% there, when omitted. If also N is omitted the original 100 dimensional
% realization is created.
%
% The code is a sparse adaption of the original dense generator function by
% Serkan Gugercin.
%

%
% MIT License
%
% Copyright (c) 2026 The MORwiki Community
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%

if nargin < 1
    N = 100;
end

if nargin < 2
    % parameter values from DOI:10.1080/00207170410001713448
    R    = 0.1;
    Rbar = 1;
    Cap  = 0.1;
    L    = 0.1;
end

z1 =  1 / L                         * ones(N, 1);
z2 = -1 / L * Rbar * R / (Rbar + R) * ones(N, 1);
z3 = -1 / L * Rbar / (Rbar + R)     * ones(N, 1);

% Every second row needs to be adapted for the second ladder step
z3(2:2:N) = -1 / Cap;
z2(1:2:N) = -1 / Cap / (Rbar + R);
z1(2:2:N) =  1 / Cap * Rbar / (Rbar + R);

% start end end are exceptional
z2(1)     = -1 / R / Cap; % sets A(1, 1)
z3(1)     = -1 / Cap;     % sets A(1, 2)
z2(N)     = -Rbar / L;    % sets A(N, N)
z1(N - 1) =  1 / L;       % sets A(N, N-1)

A = spdiags([z1 z2 z3], -1:1, N, N);
B = [(1 / R / Cap); sparse(N - 1, 1)];
C = [(-1 / R)       sparse(1, N - 1)];
D = 1 / R;

%sys_rlc = ss(A, B, C, D);