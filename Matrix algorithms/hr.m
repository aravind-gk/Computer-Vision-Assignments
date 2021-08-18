function [Q,R] = hr(A,n)
    % initialise Q and R
    Q = eye(n);
    R = A;
    
    for k = 1:n
        x = R(k:n, k);
        e = zeros(n-k+1, 1); % creating a column vector with all zeroes initially
        e(1) = 1; % first element of e contains 1, remaining all zeroes
        u = sign(x(1)) * norm(x) * e + x;
        u = u / norm(u); % normalise u
        Qk_ = eye(n-k+1) - 2 * (u * u');
        Qk = eye(n); % expand Qk_ into n-dimensions to get Qk
        Qk(k:n, k:n) = Qk_;
        Q = Q * Qk'; % find Q by multiplying it with Qk-transpose for all k in 1:n
        R(k:n, k:n) = Qk_ * R(k:n, k:n); % get R by multiplying with Qk at each step
    end
end

