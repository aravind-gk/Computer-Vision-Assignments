function [Q,R] = prop(A,n)
    % inititalising Q and R
    Q = eye(n);
    R = A;
    
    % iterating for each (i,j) in the lower triangular part of the matrix one-by-one and making it zero
    for i = 1:n
        % for each i-th column, starting from the lowest (n-th) element to (i+1)th element,
        % one-by-one make them zero using rotation matrix
        for j = n:-1:(i+1)
            % let (j-1,i)th entry 'a' and (j,i)th entry be 'b'
            a = R(j-1,i);
            b = R(j, i);
            
            % let M be the rotation matrix to perform rotation on (j-1)th and
            % j-th entries on the i-th column 
            % (thus making j-th entry equal to 0)
            cos_theta = a/sqrt(a^2 + b^2);
            sin_theta = b/sqrt(a^2 + b^2);
            M = eye(n);
            % modifying the nxn identity matrix to get desired rotating matrix
            M(j-1,j-1) = cos_theta;
            M(j-1,j) = sin_theta;
            M(j,j-1) = -sin_theta;
            M(j,j) = cos_theta;
            % performing givens rotation for each iteration
            R = M * R;
            Q = Q * M';
        end
    end
end

