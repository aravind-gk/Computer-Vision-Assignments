A = [3,3,4;3,7,6;4,6,10];
A = sym2tri(A); % making A tridiagonal

% Algorithm to convert tridiagonal A into QR form and then computing eigenvalues
n = size(A, 1);

% Find QR algorithm using Givens rotation 
% and store the maximum eigenvalue computed at the end of i-th iteration

Q = eye(n);
R = A;
maxeig = 1:(n-1); % initialising the maxeig vector
% to store the maximum eigenvalue computed at the end of i-th iteration

for i = 2:n
    % perform Givens rotation for each index (i,i-1) in A
    
    % let (i-1,i)th entry 'a' and (i,i)th entry be 'b'
    a = R(i-1,i-1);
    b = R(i,i-1);
    
    % let M be the rotation matrix to perform rotation on (i-1)th and
    % i-th entries on the (i-1)-th column 
    % (thus making i-th entry equal to 0)
    cos_theta = a/sqrt(a^2 + b^2);
    sin_theta = b/sqrt(a^2 + b^2);
    M = eye(n);
    % modifying the nxn identity matrix to get desired rotating matrix
    M(i-1,i-1) = cos_theta;
    M(i-1,i) = sin_theta;
    M(i,i-1) = -sin_theta;
    M(i,i) = cos_theta;
    % performing givens rotation for each iteration
    R = M * R;
    Q = Q * M';
    %disp(R)
    maxeig(i-1) = max(diag(R)); % storing the max eigenvalue in this iteration
end

eigen_inbuilt = max(eig(A)); % inbuilt function max eigenvalue

Xaxis = 1:(n-1);
Yaxis = 1:(n-1);
for i = 1:(n-1)
    Yaxis(i) = norm(maxeig(i) - eigen_inbuilt);
end

plot(Xaxis, Yaxis)
% Comment:
% This program for Q2_1 currently unfortunately has a few bugs and not
% giving the correct output for some test cases