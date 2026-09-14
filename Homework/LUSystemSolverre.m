function x = LUSystemSolver(A, r, b)
%LUSystemSolver  Solve Ax=b using a previously computed PA=LU factorization.
%
%   x = LUSystemSolver(A, r, b)
%
%   Assumptions (matches the common classroom LU-with-pivoting storage):
%     - A is n-by-n and stores U on and above the diagonal.
%     - The strictly lower triangle of A stores the multipliers (L entries).
%     - L is UNIT lower triangular (diagonal entries are 1 and are NOT stored).
%     - r is an (n-1)-by-1 vector recording the row swaps performed during GEPP.
%       Convention used here (and in LU_GEPP.m below): at step k, rows k and r(k)
%       were swapped. If no swap, r(k)=k.
%
%   Outputs:
%     - x solves Ax=b.
%
%   Notes:
%     - We do NOT explicitly form P, L, or U.
%     - We apply the recorded swaps directly to b to obtain Pb, then solve
%       Ly = Pb (forward substitution), and Ux = y (back substitution).

[m,n] = size(A);
if m ~= n
    error('LUSystemSolver: A must be square.');
end
if length(b) ~= n
    error('LUSystemSolver: b must have length n.');
end

% Ensure column vectors
b = b(:);
r = r(:);

% 1) Apply the same row swaps to b that were applied to A during LU_GEPP.
for k = 1:(n-1)
    pk = r(k);
    if pk ~= k
        tmp = b(k);
        b(k) = b(pk);
        b(pk) = tmp;
    end
end

% 2) Forward substitution for unit lower triangular L.
%    y(k) = b(k) - sum_{j=1}^{k-1} L(k,j)*y(j)
%    and L(k,j) is stored as A(k,j) for k>j.

y = zeros(n,1);
for k = 1:n
    s = b(k);
    for j = 1:(k-1)
        s = s - A(k,j) * y(j);
    end
    y(k) = s;  % no division because L has 1s on the diagonal
end

% 3) Back substitution for upper triangular U.
%    x(k) = ( y(k) - sum_{j=k+1}^{n} U(k,j)*x(j) ) / U(k,k)

x = zeros(n,1);
if A(n,n) == 0
    error('LUSystemSolver: zero pivot encountered in U (singular or breakdown).');
end
x(n) = y(n) / A(n,n);

for k = (n-1):-1:1
    s = y(k);
    for j = (k+1):n
        s = s - A(k,j) * x(j);
    end
    if A(k,k) == 0
        error('LUSystemSolver: zero pivot encountered in U (singular or breakdown).');
    end
    x(k) = s / A(k,k);
end

end
