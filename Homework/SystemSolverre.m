function x = SystemSolver(A, b)
%SystemSolver  Solve Ax=b using GEPP (PA=LU) followed by triangular solves.
%
%   x = SystemSolver(A, b)
%
%   Steps:
%     1) Compute PA=LU using LU_GEPP (Gaussian elimination with partial pivoting)
%     2) Solve Ly = Pb  (forward substitution)
%     3) Solve Ux = y   (back substitution)
%
%   Files used:
%     - LU_GEPP.m        (factorization)
%     - LUSystemSolver.m (solve using stored LU and swap record)

[A_lu, r] = LU_GEPPre(A);
x = LUSystemSolverre(A_lu, r, b);

end
