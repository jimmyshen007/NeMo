function n = doublefactorial(n)
%   Double FACTORIAL function.
%   DOUBLEFACTORIAL(N) for scalar even N, 
%   is the product of integers according to:

%   N!! = N*(N-2)*(N-4)...*4*2; i.e. prod(2:2:N);

%   for scalar odd N, is the product of integers according to:

%   N!! = N*(N-2)*(N-4)...*3*1; i.e. prod(1:2:N);

%   When N is an N-D matrix, DOUBLEFACTORIAL(N) is the doublefactorial for
%   each element of N.  Since double precision numbers only have about
%   15 digits, the answer is only accurate for N <= 21. For larger N,
%   the answer will have the correct order of magnitude, and is accurate for 
%   the first 15 digits.
%
%   See also PROD. and FACTORIAL
%   Date: 2006 (26 April)
%   Erick Canales Rodríguez
%  Cuban Neuroscience Center

N = n(:);
if any(fix(N) ~= N) || any(N < 0) || ~isa(N,'double') || ~isreal(N)
  error('MATLAB:doublefactorial:NNegativeInt', ...
        'N must be a matrix of non-negative integers.')
end
% n(N>170) = 171; 
% m = max([1; n(:)]);
% N = [1 1 cumprod(2:m)];
% n(:) = N(n+1);

if (ceil(n/2)==(n/2))|(n==2)  % n is even
    n = prod(2:2:n);
end

if n==1
    n=1;
elseif ceil(n/2) ~=(n/2)  % n is odd
    n = prod(1:2:n);
end
 
