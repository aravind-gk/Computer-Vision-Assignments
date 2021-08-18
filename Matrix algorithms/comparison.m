eps = 0.001; % tolerance
N = [10, 50, 100]; % matrix sizes

time1 = [0, 0, 0]; % time spent by jacobi
time2 = [0, 0, 0]; % time spent by direct method (LU partial)

for k1 = 1:3
    n = N(k1);

    % create matrix A and vector b according to question
    e = ones(n, 1);
    A = spdiags([-e 2*e -e], -1:1, n, n);
    A = full(A);
    b = rand(n, 1);
    
    tic
    % solve Ax=b using Jacobi method
    x = zeros(n, 1); % initial solution x(0)
    x_ = zeros(n, 1);
    
    % jacobi algorithm
    % we use the stopping condition: ||x_(n) - x(n)|| / ||x(n)|| < eps
    % let x denote x(k) and x_ denote x(k+1)
    iter = 0;
    
    while true
        for i = 1:n
            x_(i) = 0;
            for j = 1:n
                if j ~= i
                    x_(i) = x_(i) + A(i,j)*x(j);
                end
            end
            x_(i) = x_(i) / (-A(i,i));
            x_(i) = x_(i) + b(i)/A(i,i);
        end
        iter = iter + 1;
        
        if norm(x_ - x)/norm(x_) < eps
            break
        end
        
        x = x_;
    end
    time1(k1) = toc;
 
    % For me, the best performing direct method was LU decomposition with pivoting
    %it's analysis is as follows
    
    tic 
    
    % first find the permutation matrix P
    P = eye(n);

    % keep swapping rows in I and A such that the largest element of 
    % each column of M is placed on the topmost position of each submatrix
    swaps = 0; % variable to keep track of the number of swaps required for P
    % for each column i, do the following
    for i = 1:n
        max_row = i;
        % find the row which contains largest value in the submatrix
        % j denotes the row number 
        for j = i:n
            if A(j,i) > A(max_row, i)
                max_row = j;
            end
            if max_row ~= i
                % swap is required
                X = [i, max_row]; % rows to be swapped
                P(X,:) = P(X([2,1]),:); % syntax to swap two rows in matlab
                swaps = swaps + 1;
            end
        end
    end
    
    % Now find L and U
    % let M = PA
    M = P * A;
    L = eye(n);
    U = eye(n);
    tol = 0.0000001; % tolerance value (values smaller than this are considered zero)
    
    % Check whether LU exists and compute L and U if they exist
    LU_exists = true;
    for i=1:n
        for j = 1:n
            if i > j
                % compute L(i,j)
                L(i, j) = M(i, j);
                for k = 1:(j-1)
                    L(i, j) = L(i, j) - (L(i, k) * U(k, j));
                end
                if abs(U(j, j)) < tol
                    % division by zero, hence LU does not exist
                    LU_exists = false;
                    %disp('LU does not exist')
                    break
                else
                    L(i, j) = L(i, j) / U(j, j);
                end
            else
                % compute U(i, j)
                U(i, j) = M(i, j);
                for k = 1:(i - 1)
                    U(i, j) = U(i, j) - (L(i, k)*U(k, j));
                end
            end
        end
        if LU_exists == false
            break
        end
    end
    
    % now solve for LU = Pb
    b = P * b;
    
    % First solve Ly = b, taking Ux as y
    y = b; % initialise y as a column vector b with n elements
    
    for i = 1:n
        y(i) = b(i);
        for j = 1:(i-1)
            y(i) = y(i) - L(i,j)*y(j);
        end
    end
    
    % now solve for Ux = y
    x = y; % initialise x as a column vector of n elements
    
    for i = n:-1:1
        x(i) = y(i);
        for j = (i+1):n
            x(i) = x(i) - U(i,j)*x(j);
        end
        x(i) = x(i) / U(i, i);
    end

    time2(k1) = toc;
    
end
% Performance analysis and remarks
disp('Time using jacobi for = 10,50,100')
disp(time1)
disp('Time using LU partial for = 10,50,100')
disp(time2)

disp('Conclusion:')
disp('Even though Jacobi performs well for small values like n = 10, still the direct methods like LU partial manage to outperform Jacobi method easily for larger n like 100')
disp('Hence for very large matrices and good precision requirements, I can conclude that direct methods like LU tend to perform better than iterative methods like Jacobi')