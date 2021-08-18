%% Q 2.3.1 (a)
estimated = zeros(2); %estimated[i] contains the analytical estimate of iterations for problem i for Jacobi and GS
actual = zeros(2); %actual[i] contains the actual number of iterations for problem i for Jacobi and GS

%% first system of linear equations given in question
A = [3,1,1;1,7,3;2,0,4]; 
b = [5,11,5]';
n = 3; % n = no. of rows in square matrix A

% compute necessary matrices
% initialising L, D and U
L = zeros(n);
D = zeros(n);
U = zeros(n);

% computing L (lower triangular matrix)
for i = 1:n
    for j = 1:(i-1)
        L(i,j) = A(i,j);
    end
end

% computing U (upper triangular matrix)
for i = 1:n
    for j = (i+1):n
        U(i,j) = A(i,j); 
    end
end

% compute D (diagonal matrix)
for i = 1:n
    D(i,i) = A(i,i);
end

% compute P1 = -D-1 *(L+U) for Jacobi method
% not sure if we are allowed to use INVERSE function here, but here goes nothing
P1 = -inv(D) * (L + U);

% compute P2 = -(L+D)-1*U for Gauss siedel
P2 = -inv(L+D) * U;

% compute largest singular value of P
[~, S1, ~] = svd(P1);
sigma1 = S1(1,1);
[~, S2, ~] = svd(P2);
sigma2 = S2(1,1);

% eps is the maximum relative error (tolerance)
eps = 0.0001;

% see whether A can be solved by Jacobi or Gauss siedel methods:
if sigma1 < 1
    disp('x(k) converges to x* using Jacobi for all initial x(0)')
    
    % find expected number of iterations to converge
    k1 = ceil(log(eps) / log(sigma1));
    disp('expected number of iterations:')
    disp(k1)
else
    disp('x(k) may/may not converge to x* using Jacobi for all initial x(0)')
end

if sigma2 < 1
    disp('x(k) converges to x* using Gauss siedel for all initial x(0)')
    % find expected number of iterations to converge
    k2 = ceil(log(eps) / log(sigma2));
    disp('expected number of iterations:')
    disp(k2)
else
    disp('x(k) may/may not converge to x* using gauss siedel for all initial x(0)')
end

% storing these values for Q 2.3.1 (c)
estimated(1,1) = k1; 
estimated(1,2) = k2;
%% second system of linear equations

A = [1,5,1;9,3,3;2,1,4]; 
b = [7,15,7]';
n = 3; % n = no. of rows in square matrix A

% compute necessary matrices
% initialising L, D and U
L = zeros(n);
D = zeros(n);
U = zeros(n);

% computing L (lower triangular matrix)
for i = 1:n
    for j = 1:(i-1)
        L(i,j) = A(i,j);
    end
end

% computing U (upper triangular matrix)
for i = 1:n
    for j = (i+1):n
        U(i,j) = A(i,j); 
    end
end

% compute D (diagonal matrix)
for i = 1:n
    D(i,i) = A(i,i);
end

% compute P1 = -D-1 *(L+U) for Jacobi method
% not sure if we are allowed to use INVERSE function here, but here goes nothing
P1 = -inv(D) * (L + U);

% compute P2 = -(L+D)-1*U for Gauss siedel
P2 = -inv(L+D) * U;

% compute largest singular value of P
[~, S1, ~] = svd(P1);
sigma1 = S1(1,1);
[~, S2, ~] = svd(P2);
sigma2 = S2(1,1);

% eps is the maximum relative error (tolerance)
eps = 0.0001;

% see whether A can be solved by Jacobi or Gauss siedel methods:
if sigma1 < 1
    disp('x(k) converges to x* using Jacobi for all initial x(0)')
    
    % find expected number of iterations to converge
    k1 = ceil(log(eps) / log(sigma1));
    disp('expected number of iterations:')
    disp(k1)
else
    disp('x(k) may/may not converge to x* using Jacobi for all initial x(0)')
end

if sigma2 < 1
    disp('x(k) converges to x* using Gauss siedel for all initial x(0)')
    % find expected number of iterations to converge
    k2 = ceil(log(eps) / log(sigma2));
    disp('expected number of iterations:')
    disp(k2)
else
    disp('x(k) may/may not converge to x* using gauss siedel for all initial x(0)')
end

%% second system of equations (after row permutation)
disp('Since iterative methods for the second system of equations is not converging, we wil permute the rows of A to make it converge')
A = [9,3,3;1,5,1;2,1,4]; % row 1 and row 2 exchanged
b = [15,7,7]';
n = 3; % n = no. of rows in square matrix A

% compute necessary matrices
% initialising L, D and U
L = zeros(n);
D = zeros(n);
U = zeros(n);

% computing L (lower triangular matrix)
for i = 1:n
    for j = 1:(i-1)
        L(i,j) = A(i,j);
    end
end

% computing U (upper triangular matrix)
for i = 1:n
    for j = (i+1):n
        U(i,j) = A(i,j); 
    end
end

% compute D (diagonal matrix)
for i = 1:n
    D(i,i) = A(i,i);
end

% compute P1 = -D-1 *(L+U) for Jacobi method
% not sure if we are allowed to use INVERSE function here, but here goes nothing
P1 = -inv(D) * (L + U);

% compute P2 = -(L+D)-1*U for Gauss siedel
P2 = -inv(L+D) * U;

% compute largest singular value of P
[~, S1, ~] = svd(P1);
sigma1 = S1(1,1);
[~, S2, ~] = svd(P2);
sigma2 = S2(1,1);

% eps is the maximum relative error (tolerance)
eps = 0.0001;

% see whether A can be solved by Jacobi or Gauss siedel methods:
if sigma1 < 1
    disp('x(k) converges to x* using Jacobi for all initial x(0)')
    
    % find expected number of iterations to converge
    k1 = ceil(log(eps) / log(sigma1));
    disp('expected number of iterations:')
    disp(k1)
else
    disp('x(k) may/may not converge to x* using Jacobi for all initial x(0)')
end

if sigma2 < 1
    disp('x(k) converges to x* using Gauss siedel for all initial x(0)')
    % find expected number of iterations to converge
    k2 = ceil(log(eps) / log(sigma2));
    disp('expected number of iterations:')
    disp(k2)
else
    disp('x(k) may/may not converge to x* using gauss siedel for all initial x(0)')
end

% storing these values for Q 2.3.1 (c)
estimated(2,1) = k1; 
estimated(2,2) = k2;

%% Q 2.3.1 (b) (Problem 1)

% These N and y's will be used in part (c) of the question, and will be explained later
N = 1:5;
y11 = [0, 0, 0, 0, 0];
y12 = [0, 0, 0, 0, 0];
y21 = [0, 0, 0, 0, 0];
y22 = [0, 0, 0, 0, 0];

% Enter any matrix in place of A and b, and update n accordingly
A = [3,1,1;1,7,3;2,0,4]; 
b = [5,11,5]';
n = 3; % n = no. of rows in square matrix A

x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);
eps = 0.001;

% jacobi algorithm
% we use the stopping condition: ||Ax - b|| < eps
% let x denote x(k) and x_ denote x(k+1)
iter = 0;

while norm(A*x - b) > eps
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
    if iter <= 5
        y11(iter) = norm(x_ - x)/norm(x_); 
    end
    x = x_;
    if iter > 1000
        % to prevent infinite loops
        disp('not converging')
        break
    end
end

disp('solution of Ax=b (problem 1) using Jacobi method is:')
disp(x)
disp('#of iterations using Jacobi:')
disp(iter)

% storing these values for Q 2.3.1 (c)
actual(1,1) = iter;
%% gauss-siedel algorithm

% we use the stopping condition: ||Ax - b|| < eps
% let x denote x(k) and x_ denote x(k+1)
iter = 0;
x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);

while norm(A*x - b) > eps
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
    if iter <= 5
        y12(iter) = norm(x_ - x)/norm(x_); 
    end
    x = x_;
    if iter > 1000
        % to prevent infinite loops
        disp('not converging')
        break
    end
end

disp('solution of Ax=b (problem 1) using gauss siedel method is:')
disp(x)
disp('#of iterations using gauss siedel :')
disp(iter)

% storing these values for Q 2.3.1 (c)
actual(1,2) = iter;
%% Q 2.3.1 (b) (Problem 2)

% Enter any matrix in place of A and b, and update n accordingly
A = [9,3,3;1,5,1;2,1,4]; 
b = [15,7,7]';
n = 3; % n = no. of rows in square matrix A

x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);
eps = 0.001;

% jacobi algorithm
% we use the stopping condition: ||Ax - b|| < eps
% let x denote x(k) and x_ denote x(k+1)
iter = 0;

while norm(A*x - b) > eps
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
    if iter <= 5
        y21(iter) = norm(x_ - x)/norm(x_); 
    end
    x = x_;
    if iter > 1000
        % to prevent infinite loops
        disp('not converging')
        break
    end
end

disp('solution of Ax=b (problem 2) using Jacobi method is:')
disp(x)
disp('#of iterations using Jacobi:')
disp(iter)

% storing these values for Q 2.3.1 (c)
actual(2,1) = iter;
%% gauss-siedel algorithm

% we use the stopping condition: ||Ax - b|| < eps
% let x denote x(k) and x_ denote x(k+1)
iter = 0;
x = zeros(n, 1); % initial solution x(0)
x_ = zeros(n, 1);

while norm(A*x - b) > eps
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
    %disp(norm(x_ - x)/norm(x_))
    if iter <= 5
        y22(iter) = norm(x_ - x)/norm(x_); 
    end
    x = x_;
    if iter > 1000
        % to prevent infinite loops
        disp('not converging')
        break
    end
end


disp('solution of Ax=b (problem 2) using gauss siedel method is:')
disp(x)
disp('#of iterations using gauss siedel :')
disp(iter)

% storing these values for Q 2.3.1 (c)
actual(2,2) = iter;
%% Q 2.3.1 (c)
fprintf('\nestimated iterations for problem 1 using jacobi = %d\n', estimated(1,1))
fprintf('estimated iterations for problem 1 using gauss siedel = %d\n', estimated(1,2))
fprintf('estimated iterations for problem 2 using jacobi = %d\n', estimated(2,1))
fprintf('estimated iterations for problem 2 using gauss siedel = %d\n', estimated(2,2))

fprintf('\nactual iterations for problem 1 using jacobi = %d\n', actual(1,1))
fprintf('actual iterations for problem 1 using gauss siedel = %d\n', actual(1,2))
fprintf('actual iterations for problem 2 using jacobi = %d\n', actual(2,1))
fprintf('actual iterations for problem 2 using gauss siedel = %d\n', actual(2,2))

fprintf('\nAccording to estimated values, Gauss-siedel converges sligtly faster as compared to Jacobi method.\n')
fprintf('However, upon actual execution of both the algorithms, Gauss-siedel algorithm converges much quicker, almost twice as fast as Jacobi algorithm.\n')

% Plot to demonstrate convergence of jacobi and gauss siedel
% y11: problem 1, jacobi
% y12: problem 1, gauss
% y21: problem 2, jacobi
% y22: problem 2, gauss

%% plotting graphs to show convergence (relative error) in problem 1
% (here, the blue plot denotes jacobi and orange denotes gauss

%disp(y11)
%disp(y12)
%disp(y21)
%disp(y22)
%disp(N)

plot(N, y11)
title('Problem 1 convergence (relative error)')
hold on
plot(N, y12)
hold off

%% plotting graphs to show convergence (relative error) in problem 2
% (here, the blue plot denotes jacobi and orange denotes gauss

plot(N, y21)
title('Problem 2 convergence (relative error)')
hold on
plot(N, y22)
hold off