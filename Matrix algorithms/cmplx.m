%% Either Take matrix A as user input 
n=input('order of your square matrix is? ');
for i=1:n^2
    A(i)=input('elements-');
end
A=reshape(A,n,n)';

%% Or Take matrix A as hardcoded input
n = 3;
A = [1, 2, 3; 4, 5, 6; 7, 8, 9];
%A = zeros(3, 3); %(case where LU does not exist)
A=reshape(A,n,n); 

% take 'b' as input
b = [5, 6, 7]';
%% Create L and U matrices (initialise)
L = eye(n);
U = eye(n);
tol = 0.0000001; % tolerance value (values smaller than this are considered zero)

%% Find the LU decomposition of A (without pivot)
% Check whether LU exists and compute L and U if they exist
LU_exists = true;
for i=1:n
    for j = 1:n
        if i > j
            % compute L(i,j)
            L(i, j) = A(i, j);
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
            U(i, j) = A(i, j);
            for k = 1:(i - 1)
                U(i, j) = U(i, j) - (L(i, k)*U(k, j));
            end
        end
    end
    if LU_exists == false
        break
    end
end

%% Display L and U

if LU_exists == true
    disp(L) 
    disp(U)
else
    disp('LU decomposition does not exist')
end

for i = 1:n
    if abs(U(i,i)) < tol
        % condition equivalent to saying U(i,i) == 0
        disp('The diagonal of U contains zeros because the matrix A is singular.')
        break
    end
end



%% Verify that LU decomposition has been performed correctly
if LU_exists == true
    disp(L*U)
else
    disp('LU decomposition does not exist')
end

%% Solve Ax = b using LU decomposition, equivalently solve LUx = b

% (solve Ax=b only if LU actually exists)
if LU_exists == true
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
        if abs(U(i,i)) < tol
            % condition equivalent to saying "if U(i,i) == 0"
            % infinite solutions possible in this case
            % so pick any suitable value of x
            x(i) = 1;
        else
            x(i) = y(i);
            for j = (i+1):n
                x(i) = x(i) - U(i,j)*x(j);
            end
            x(i) = x(i) / U(i, i);
        end
    end
else
    disp('LU decomposition does not exist')
end

%% display the obtained solution x
if LU_exists == true
    disp(x)
else
    disp('LU decomposition does not exist')
end

%% Compare Ax with b, to verify that Ax = b 
if LU_exists == true
    disp([A*x, b])
else
    disp('LU decomposition does not exist')
end 

%% Verify that LUx = b
if LU_exists == true
    disp([L*U*x, b])
else
    disp('LU decomposition does not exist')
end