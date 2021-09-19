%% input matrix A 

A = [1,2,3;4,5,6;7,8,9]; % can change the value of matrix A here, and update 'n' as well
n = 3;

%% first find QR decomposition of A using householder reflections

[Q, R] = hr(A, n);
fprintf("Q:")
disp(Q)
fprintf("R:")
disp(R)

%% get LU decomposition from QR
% basically, find the LU decomposition of Q to get L and U
% and let U1 = U * R, and product of 2 upper triangular matrices is upper triangular
% so L and U1 is a valid LU decomposition of A = QR

L = eye(n);
U = eye(n);
tol = 0.0000001; % tolerance value (values smaller than this are considered zero)

for i=1:n
    for j = 1:n
        if i > j
            % compute L(i,j)
            L(i, j) = Q(i, j);
            for k = 1:(j-1)
                L(i, j) = L(i, j) - (L(i, k) * U(k, j));
            end
            L(i, j) = L(i, j) / U(j, j);
        else
            % compute U(i, j)
            U(i, j) = Q(i, j);
            for k = 1:(i - 1)
                U(i, j) = U(i, j) - (L(i, k)*U(k, j));
            end
        end
    end
end

%% LU decomposition of Q
disp(L)
disp(U)
disp(L*U)

%% Get LU decomposition of A
U1 = U * R;

fprintf("L:")
disp(L)
fprintf("U1:")
disp(U1)

% verify that L * U1 = A
fprintf("L * U1:")
disp(L * U1)
