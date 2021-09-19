%% case when n = 2 
mu = 0; % mean of the required normal distribution, 0 in this case
C = [0.025, 0.0075; 0.0075, 0.007]; % The covariance matrix 
L = cholprog(C); % cholesky decomposition of covariance matrix 

fprintf('Cholesky decomposition of covariance matrix:')
disp(L)

% generate correlated normal random vectors using L

M = 1000; % number of desired random vectors
Y = corrNRV(L, M, mu);

Y1 = Y(1,:);
Y2 = Y(2,:);
scatter(Y1, Y2)

%% case when n = 3
mu = 0; % mean of the required normal distribution, 0 in this case 
C = [0.025, 0.0075, 0.00175; 0.0075, 0.007, 0.00135; 0.00175, 0.00135, 0.00043]; % The covariance matrix 
L = cholprog(C); % cholesky decomposition of covariance matrix 

fprintf('Cholesky decomposition of covariance matrix:')
disp(L)

% generate correlated normal random vectors using L

M = 1000; % number of desired random vectors
Y = corrNRV(L, M, mu);

Y1 = Y(1,:);
Y2 = Y(2,:);
Y3 = Y(3,:);
scatter3(Y1, Y2, Y3)