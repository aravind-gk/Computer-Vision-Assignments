%% initialise the matrices
n = 10; % dimension of the square matrices

% generate two random vectors
v1 = randn(n, 1);
v2 = randn(n, 1);

% make v1 and v2 orthogonal to each other, then make them orthonormal
v1 = v1 - ((v1'*v2)/(v2'*v2))*v2;
v1 = v1/norm(v1);
v2 = v2/norm(v2);

% Generate A using v1 and v2
A = 50000*(v1*v1') + 2*(v2*v2');

% generate x and b
x = randn(n, 1);
b = A*x;

% let x be the actual x and x_ be the estimated x
%% now solve Ax_ = b, or equivalently, Rx_ = Q'b

% first find the QR decomposition of A using the three methods
[Q1, R1] = gs(A, n);
[Q2, R2] = hr(A, n);
[Q3, R3] = prop(A, n);

% Q_ and R_ are arrays of Q and R, using each of the 3 methods
% Now we will solve for each Q and R and record the error
Q_ = {Q1, Q2, Q3}; 
R_ = {R1, R2, R3};

error = [];

for k = 1:3
    Q = Q_{k};
    R = R_{k};
    
    % Now solve Rx_ = Q'b
    % let c = Q'b, then solve for Rx_ = c using back-substitution of upper triangular matrix R
    c = Q' * b;
    x_ = zeros(n, 1);
    
    for i = n:-1:1
        x_(i) = c(i);
        for j = (i+1):n
            x_(i) = x_(i) - R(i,j)*x_(j);
        end
        x_(i) = x_(i)/R(i,i);
    end
    
    % find the error between x and x_ and then append it to the "error" vector
    error = [error, norm(x - x_)];
end

%% display the errors by all 3 methods

method = {'Gram schmidt ', 'householder reflections ', 'rotational (prop) '};
disp('the errors for each method are:')
for i = 1:3
    disp(method{i})
    
    disp(error(i))
end

% find the best method (with least error)
[~, index] = min(error); 
disp('The best method is: ')
disp(method{index}) 

%% temp
[Q1, R1] = gs(A, n);
[Q2, R2] = hr(A, n);
[Q3, R3] = prop(A, n);
Q = {Q1, Q2, Q3};

%% temp2
disp(Q{2})
