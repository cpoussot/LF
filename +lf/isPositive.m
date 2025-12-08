% <strong>MDSPACK</strong> (A > 0)
% MDSPACK.TOOLS.ISPOSITIVE - Check matrix positive definitness 
%
% Syntax
%  [tag, eigA] = mdspack.tools.isPositive(A)
%  [tag, eigA] = mdspack.tools.isPositive(A, tol)
%  
% Input arguments 
%  - A : square matrix 
% 
% Output arguments
%  - tag  : positivity tag (boolean)
%  - eigA : eigenvalues of A
% 
% Optional arguments
%  - tol : tolerance of the minimal real part of eig(A) (positive real, 
%          default is 1e-10)
% 
% Description
% Checks if A > 0. Numerically it checks min(real(eig(A))) > -tol.
% 

function [t,eigA] = isPositive(A,tol)

if nargin < 2
    tol = 1e-10;
end
eigA    = eig(A);
t       = false;
if min(real(eigA)) > -tol
    t = true;
end
