%% defining the matrices and taking input

n = 10; % n is the size of the square matrix

% define the first col of the matrix
col = zeros(n, 1);
col(1) = 2;
col(2) = -1;

% define b
b = zeros(n, 1);
b(1) = 11;

%% construct toplitz matrix A from col

A = zeros(n);

for i = 1:n
    A(i,1) = col(i);
    A(1,i) = col(i);
end

for i = 2:n
    for j = 2:n
        A(i,j) = A(i-1,j-1);
    end
end

disp('the toplitz matrix A (10x10) is as follows:')
disp(A)

%% Solve Ax=b using Jacobi method

x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);
eps = 0.0001;

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

disp('solution of Ax=b using Jacobi method is:')
disp(x)
disp('#of iterations using Jacobi:')
disp(iter)

%% Solve Ax=b using gauss-siedel

% let x denote x(k) and x_ denote x(k+1)
iter = 0;
x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);
eps = 0.0001;

while true
    for i = 1:n
        x_(i) = 0;
        for j = 1:(i-1)
            x_(i) = x_(i) + A(i,j)*x_(j);
        end
        for j = (i+1):n
            x_(i) = x_(i) + A(i,j)*x(j);
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

disp('solution of Ax=b using gauss siedel method is:')
disp(x)
disp('#of iterations using gauss siedel :')
disp(iter)

%% repeat the above problem using 20*20 toplitz matrix 

n = 20; % n is the size of the square matrix

% define the first col of the matrix
col = zeros(n, 1);
col(1) = 2;
col(2) = -1;

% define b
b = zeros(n, 1);
b(1) = 21;

%% construct toplitz matrix A from col

A = zeros(n);

for i = 1:n
    A(i,1) = col(i);
    A(1,i) = col(i);
end

for i = 2:n
    for j = 2:n
        A(i,j) = A(i-1,j-1);
    end
end

disp('the toplitz matrix A (20x20) is as follows:')
disp(A)

%% Solve Ax=b using Jacobi method

x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);
eps = 0.0001;

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

disp('solution of Ax=b using Jacobi method is:')
disp(x)
disp('#of iterations using Jacobi:')
disp(iter)

%% Solve Ax=b using gauss-siedel

% let x denote x(k) and x_ denote x(k+1)
iter = 0;
x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);
eps = 0.0001;

while true
    for i = 1:n
        x_(i) = 0;
        for j = 1:(i-1)
            x_(i) = x_(i) + A(i,j)*x_(j);
        end
        for j = (i+1):n
            x_(i) = x_(i) + A(i,j)*x(j);
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

disp('solution of Ax=b using gauss siedel method is:')
disp(x)
disp('#of iterations using gauss siedel :')
disp(iter)