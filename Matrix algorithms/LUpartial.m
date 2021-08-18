%% Take matrix A as hardcoded input
% (You may change these values and test for other input matrices)

n = 3;
A = [1, 2, 3; 4, 5, 6; 7, 8, 9];

% take 'b' as input
b = [5, 6, 7]';
disp(A)

%% find the permutation matrix P such that LU decomposition exists for PA

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

disp('The permutation matrix P is as follows:')
disp(P)
disp('Now finding the LU decomposition of PA')

%% LU decomposition (code copied from cmplx.m)
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

%% display L and U
fprintf("L:")
disp(L)
fprintf("U:")
disp(U)

%% verify that PA = LU, or, A = P'LU
disp(P' * L * U)

%% determinant of A = P'LU

detP = (-1)^swaps; % determinant of permutation matrix
detL = 1; % since all diagonal entries are 1

detU = 1;
for i = 1:n
    detU = detU * U(i,i);
end

fprintf('determinant of A = P"LU is:')
disp(detP * detL * detU)

%% inverse of A = P'LU
%{
% let the inverse of P*A = L*U be A1
A1 = zeros(n);
I = eye(n);

for k = 1:n
    % solve for each i-th column
    % first solve Ly 
    % First solve Ly = b, taking Ux as y
    b = I(1:n,k); % b is the ith column of identity matrix
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
    
    % now x is the k-th column of the inverse of LU
    A(1:n,k) = x;
end

% to get inverse of A, we do A1 = A1 * P
% now A1 is the inv of A
A1 = A1 * P;
%}
