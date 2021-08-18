T = [0.1, 0.01, 0.001, 0.0001, 0.00001]; % tolerance levels
N = [10, 50, 100]; % size of matrix 

it = zeros(3, 5);

% for every pair of tolerance and n, run the algorithm and record the no.
% of iterations
for k1 = 1:3
    for k2 = 1:5
        n = N(k1);
        eps = T(k2);
        
        % create matrix A and vector b according to question
        e = ones(n, 1);
        A = spdiags([-e 2*e -e], -1:1, n, n);
        A = full(A);
        b = rand(n, 1);
        
        % solve Ax=b using Jacobi method
        x = zeros(n, 1); % initial solution x(0)
        x_ = zeros(n, 1);
        
        % jacobi algorithm
        % we use the stopping condition: ||x_(n) - x(n)|| / ||x(n)|| < eps
        % let x denote x(k) and x_ denote x(k+1)
        iter = 0;
        
        while true
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
            
            if norm(x_ - x)/norm(x_) < eps
                break
            end
            
            x = x_;
        end
        
        it(k1,k2) = iter;
    end
end

%% 
disp('the table of list of iterations for different values of n and tolerance:')
disp('Rows are different values of n and columns are different tolerance')
disp(it)

%% plot the results
% using a loglog plot as normal plots are too skewed
disp('blue plot denotes n = 10')
disp('red plot denotes n = 50')
disp('orange plot denotes n = 100')

for i = 1:3
    y = it(i,:);
    loglog(T, y)
    if i == 1
        hold on
    end
end
hold off

%% 

disp('Comments:')
disp('The number of iterations required increases exponentially with decrease in tolerance')
disp('The iterations also increases rapidly with increase in n')
disp('Hence jacobi method is not very suitable for very large and low tolerance systems')

