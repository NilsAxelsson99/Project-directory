function [V, D] = compute_eigen(A)
%COMPUTE_EIGEN Computes eigenvalues and eigenvectors of a matrix
%
% Inputs:
%   A - Square matrix (n x n)
%
% Outputs:
%   V - Matrix whose columns are eigenvectors
%   D - Diagonal matrix of eigenvalues

    % Check that the matrix is square
    [m,n] = size(A);
    if m ~= n
        error('Input matrix must be square.');
    end

    % Compute eigenvalues and eigenvectors
    [V,D] = eig(A);

end