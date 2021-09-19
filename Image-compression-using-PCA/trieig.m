function [eigen] = trieig(M)
%TRIEIG Summary of this function goes here
    A = M; % making a copy of the matrix M
    % assuming that input M is already a tridiagonal matrix
    
    % Algorithm to convert tridiagonal A into QR form and then computing eigenvalues
    n = size(A, 1);
    
    % Find QR algorithm using Givens rotation 
    % Givens rotation is performed only on the first subdiagonal elements 
    % thus only (n-1) rotations are needed to be performed for tridiagonal A
    
    Q = eye(n);
    R = A;

    for i = 2:n
        % perform Givens rotation for each index (i,i-1) in A
        
        % let (i-1,i)th entry 'a' and (i,i)th entry be 'b'
        a = R(i-1,i-1);
        b = R(i,i-1);
        
        % let M be the rotation matrix to perform rotation on (i-1)th and
        % i-th entries on the (i-1)-th column 
        % (thus making i-th entry equal to 0)
        cos_theta = a/sqrt(a^2 + b^2);
        sin_theta = b/sqrt(a^2 + b^2);
        M = eye(n);
        % modifying the nxn identity matrix to get desired rotating matrix
        M(i-1,i-1) = cos_theta;
        M(i-1,i) = sin_theta;
        M(i,i-1) = -sin_theta;
        M(i,i) = cos_theta;
        % performing givens rotation for each iteration
        R = M * R;
        Q = Q * M';
        %disp(R)
        % NOTE: Since for each iteration, a matrix multiplication
        % (rotation) is performed, the per-iteration cost of this algorithm
        % is O(n^3) 
        % although it could be a bit more optimized since most of the
        % entries in this matrix are zero
    end
    fprintf('R from trieig function:')
    disp(R)
    eigen = diag(R); % the eigenvalues of A=QR are the diagonal elements of R
end