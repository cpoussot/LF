function Rxy = tcorrelation(x, y, key)

if nargin < 2
    y = x;
end
if nargin < 3
    key = 'B';
end
% x = x - mean(x);
% y = y - mean(y);
N = length(x);

switch key
    case 'B' % biased 
        alpha = @(m) 1/N;
    case 'U' % un-biased
        alpha = @(m) 1/(N-m);
end
%
Rxy = zeros(N,1);
for m = 1:N
    Rxy(m) = alpha(m-1) * x(1,1:N-m)*y(1,m+1:N)';
    %Rxy(m) = alpha(m-1) * x(1,1:N-m)'*y(1,m:N);
end
