%% Enter the value of 'n' (dimension of square matrix) (by default, 3)
n = 71;

% Initialise matrix A and Identity matrix I
A = 0.00001 * eye(n) + hilb(n);
I = eye(n);

%% Apply QR decomposition using each of the 3 algorithms and note the error between QtQ and I
[Q, ~] = gs(A, n);
error_gs = frob(Q'*Q, I, n);

[Q, ~] = hr(A, n);
error_hr = frob(Q'*Q, I, n);  

[Q, ~] = prop(A, n);
error_prop = frob(Q'*Q, I, n);

% display all errors
disp('errors between Qt*Q and I are:')
disp('error gs:')
disp(error_gs)
disp('error hr:')
disp(error_hr)
disp('error prop:')
disp(error_prop)
%% find the algorithm which gives least error between QR and A

[Q, R] = gs(A, n);
error_gs = frob(Q*R, A, n);

[Q, R] = hr(A, n);
error_hr = frob(Q*R, A, n);  

[Q, R] = prop(A, n);
error_prop = frob(Q*R, A, n);

% display all errors
disp('errors between Q*R and A are:')
disp('error gs:')
disp(error_gs)
disp('error hr:')
disp(error_hr)
disp('error prop:')
disp(error_prop)

% find and display min-error 
error = [error_gs, error_hr, error_prop];
[min_error, index] = min(error);

if index == 1
    best_QR = 'Gram-Schmidt!';
elseif index == 2
    best_QR = 'Householder reflections!';
else 
    best_QR = 'rotation transformation!';
end

disp('The best QR decomposition algorithm is ')
disp(best_QR)
disp('and the least error (between QR and A) achieved in this case is ')
disp(min_error)

%% display all errors
disp('error gs:')
disp(error_gs)
disp('error hr:')
disp(error_hr)
disp('error prop:')
disp(error_prop)

%% local function to compute frobenius norm
function output = frob(A, B, n)
    output = 0;
    for i = 1:n
        for j = 1:n
            output = output + (A(i,j) - B(i,j))^2;
        end
    end
    output = sqrt(output);
end