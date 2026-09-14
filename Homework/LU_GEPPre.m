function [A, r] = LU_GEPP(A)
%LU_GEPP  Compute PA=LU using Gaussian Elimination with Partial Pivoting.
%
%   [A, r] = LU_GEPP(A)
%
%   On output, A stores both L and U in-place:
%     - U is stored on and above the diagonal of A
%     - L multipliers are stored below the diagonal of A
%     - L is unit lower triangular (diagonal = 1, not stored)
%
%   r is an (n-1)-by-1 vector recording the row swaps:
%     at step k, row k was swapped with row r(k).
%     If no swap occurred, r(k) = k.
%
%   This format is designed to pair with LUSystemSolver.m.

[m,n] = size(A);
if m ~= n
    error('LU_GEPP: A must be square.');
end

r = (1:n-1)';  % will be overwritten with actual pivot row choices

for k = 1:(n-1)

    % --- Partial pivoting: find pivot row among k..n that maximizes |A(i,k)|
    pivotRow = k;
    pivotVal = abs(A(k,k));
    for i = (k+1):n
        if abs(A(i,k)) > pivotVal
            pivotVal = abs(A(i,k));
            pivotRow = i;
        end
    end

    % Record the pivot row used at step k
    r(k) = pivotRow;

    % Swap rows in A if needed
    if pivotRow ~= k
        tmp = A(k,:);
        A(k,:) = A(pivotRow,:);
        A(pivotRow,:) = tmp;
    end

    % Check pivot
    if A(k,k) == 0
        error('LU_GEPP: zero pivot encountered at step k=%d (matrix singular).', k);
    end

    % --- Elimination below the pivot
    for i = (k+1):n
        % Multiplier
        A(i,k) = A(i,k) / A(k,k);
        % Row update on the remaining entries
        for j = (k+1):n
            A(i,j) = A(i,j) - A(i,k) * A(k,j);
        end
    end

end

% Last pivot check
if A(n,n) == 0
    error('LU_GEPP: zero pivot encountered at last step (matrix singular).');
end

end
