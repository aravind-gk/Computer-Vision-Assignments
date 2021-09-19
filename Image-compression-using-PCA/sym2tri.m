function [A] = sym2tri(M)
% A is the input symmetric matrix 
% This function converts A into a tridiagonal matrix and returns it 
    A = M; % creating a copy of the matrix M 
    n = size(A, 1); % size of square matrix A
    
    % Same householder's reflections algorithm as in PGA1
    for i = 1:(n-2)
        x = A(i+1:n,i); % first row of the i-th submatrix
        e = zeros(n-i,1); % creating a column vector with all zeros
        e(1) = 1; % making the first element 1 to create the first standard basis vector
        v = sign(x(1))*norm(x)*e + x;
        v = v/norm(v); % normalising v 
        % Apply the reflection on A from both above diagonal and below diagonal
        A(i+1:n,i:n) = A(i+1:n,i:n) - 2 * v * (v' * A(i+1:n,i:n));
        A(1:n,i+1:n) = A(1:n,i+1:n) - 2 * (A(1:n,i+1:n) * v) * v';
    end 
    % return the hessenberg tridiagonal matrix A
end