function [Q,R] = gs(A,n)
    Q = zeros(n);
    R = zeros(n);
    
    % u is a temporary matrix to store intermediate values
    u = zeros(n);
    
    for i = 1:n
        a = A(:,i);
        u(:,i) = a;
        % subtract the projections from other vectors to create an orthogonal vector
        for j = 1:(i-1)
            u(:,i) = u(:,i) - ((u(:,j)'*a)/(u(:,j)'*u(:,j))) * u(:,j);
        end
        % normalise to create orthonormal vector
        Q(:,i) = u(:,i)/norm(u(:,i));
    end
    
    % since Q is an orthogonal matrix, A = QR implies Q'A = R
    % use this fact to find R
    for i = 1:n
        for j = i:n
            R(i,j) = Q(:,i)'*A(:,j);
        end
    end
end

