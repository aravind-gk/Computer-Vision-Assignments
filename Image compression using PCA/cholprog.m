function [L] = cholprog(M)
%CHOLPROG Summary of this function goes here
    A = M; % make a copy of the matrix M
    n = size(A, 1);
    
    L = zeros(n); % initialise L matrix
    
    % perform cholesky decomposition such that A = L*L'
    % fill in more comments for this program
    for i = 1:n
        for j = 1:i
            sum = 0;
            if i == j
               % perform summation for diagonals
               for k = 1:(j-1)
                   sum = sum + L(j,k)*L(j,k);
               end
               L(j,j) = sqrt(A(j,j) - sum);
            else
                for k = 1:(j-1)
                    sum = sum + (L(i,k) * L(j,k));
                end
                L(i,j) = (A(i,j) - sum)/L(j,j);
            end
        end
    end
end

