% <strong>MDSPACK</strong> (A = A^H)
% MDSPACK.TOOLS.ISHERMITIAN - Check matrix is Hermitian
%
% Syntax
%  [tag, err] = mdspack.tools.isHermitian(A)
%  [tag, err] = mdspack.tools.isHermitian(A, tol)
%  
% Input arguments 
%  - A : square matrix 
% 
% Output arguments
%  - tag : Hermitian tag (boolean)
%  - err : error computed as norm(A-A^H)
% 
% Optional arguments
%  - tol : tolerance of the error (positive real, default is 1e-10)
% 
% Description
% Checks if A = A^H. Numerically it checks norm(A-A^H) < tol.
% 

function [t,h] = isHermitian(A,tol)

if nargin < 2
    tol = 1e-10;
end
h = norm(A-A');
t = false;
if h < tol
    t = true;
end