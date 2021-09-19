function [Y] = corrNRV(L, M, mu)
%CORRNRV Summary of this function goes here
% This function takes in as input the cholesky decomposition L of the covariance matrix, 
% generates 'M' uncorrelated random vectors and uses L to generate correlated random vectors 
% (which are returned as output)
% Also, 'mu' is the mean of the correlated normal random vectors

    n = size(L, 1);
    X = randn(n, M); % Here each column of the matrix 'X' is a random vector of size Mx1
    Y = L * X; % Multiplying 'X' by the cholesky decomposition of covariance matrix to get correlated random vectors of mean 0
    %Y = Y + mu; % shifting the random vectors by the given mean 'mu'
    %disp(X)
    %disp(Y)
end

