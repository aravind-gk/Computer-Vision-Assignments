function [singular] = findsing(A)
% compute the eigendecomposition of A'A and get the square roots of the eigenvalues
    M = A' * A; % compute A'A
    n = size(A, 1); % size of the matrix A
    
    % first find the hessenberg form of A'A
    H = sym2tri(M);
    
    % then compute the eigenvalues from the hessenberg
    eigen = trieig(H);
    
    % then compute the square root of the eigenvalues
    singular = 1:n;
    for i = 1:n
       singular(i) = sqrt(eigen(i)); 
    end
    
    fprintf('the singular values of given matrix are:')
    disp(singular)
end

