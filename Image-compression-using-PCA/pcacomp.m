% Program to perform Principal Component Analysis 

I = im2double(imread('watch.bmp')); % Read the image file 
% Extract the three color channels 
R = I(:,:,1); % Red Channel
G = I(:,:,2); % Green Channel
B = I(:,:,3); % Blue Channel

%% initialise blocks
XR = zeros(64, 12288); % for red
XG = zeros(64, 12288); % for green
XB = zeros(64, 12288); % for blue

% now extract and resize the blocks from 8x8 to 64x1

% for red blocks
for i = 0:95
    for j = 0:127
        x = R(8*i+1:8*(i+1), 8*j+1:8*(j+1));
        x = reshape(x, [64, 1]);
        XR(1:64,i*128+j+1) = x;
    end
end

% green blocks
for i = 0:95
    for j = 0:127
        x = G(8*i+1:8*(i+1), 8*j+1:8*(j+1));
        x = reshape(x, [64, 1]);
        XG(1:64,i*128+j+1) = x;
    end
end

% blue blocks
for i = 0:95
    for j = 0:127
        x = B(8*i+1:8*(i+1), 8*j+1:8*(j+1));
        x = reshape(x, [64, 1]);
        XB(1:64,i*128+j+1) = x;
    end
end
disp('so far so good 1')

%% calculate sample mean and covariance matrix

meanR = zeros(64,1);
meanG = zeros(64,1);
meanB = zeros(64,1);

covR = zeros(64);
covG = zeros(64);
covB = zeros(64);

for i = 1:12288
    meanR = meanR + XR(:,i);
    meanG = meanG + XG(:,i);
    meanB = meanB + XB(:,i);
end

meanR = meanR * (1/12288);
meanG = meanG * (1/12288);
meanB = meanB * (1/12288);

for i = 1:12288
    covR = covR + (XR(:,i) - meanR)*(XR(:,i) - meanR)';
    covG = covG + (XG(:,i) - meanG)*(XG(:,i) - meanG)';
    covB = covB + (XB(:,i) - meanB)*(XB(:,i) - meanB)';
end

covR = covR * (1/12288);
covG = covG * (1/12288);
covB = covB * (1/12288);

disp('so far so good 2')

%% compute eigenvalues and eigenvectors of the covariance matrices
% (the following code is from the question paper) 

% for red
[V_R, lambdaR] = eig(covR);
lambdaR = diag(lambdaR);
[lambdaR, indices] = sort(lambdaR, 'descend');
V_R = V_R(:,indices);

% for green
[V_G, lambdaG] = eig(covG);
lambdaG = diag(lambdaG);
[lambdaG, indices] = sort(lambdaG, 'descend');
V_G = V_G(:,indices);

% for blue
[V_B, lambdaB] = eig(covB);
lambdaB = diag(lambdaB);
[lambdaB, indices] = sort(lambdaB, 'descend');
V_B = V_B(:,indices);

disp('eigenvalues and eigenvectors computed and sorted')
disp('so far so good 3')

%% approximate the image using the eigenbasis by taking the 5 principal eigenvectors
K = 5;

% Let J denote the approximated image
% with JR, JG, JB denoting the RGB components respectively

% initialise with given dimensions
JR = zeros(64,12288);
JG = zeros(64,12288);
JB = zeros(64,12288);

% using the formula for constructing image using an orthobasis
for i = 1:12288
    JR(:,i) = meanR;
    JG(:,i) = meanG;
    JB(:,i) = meanB;
    
    for k = 1:K
        JR(:,i) = JR(:,i) + ((XR(:,i) - meanR)' * V_R(:,k)) * V_R(:,k);
        JG(:,i) = JG(:,i) + ((XG(:,i) - meanG)' * V_G(:,k)) * V_G(:,k);
        JB(:,i) = JB(:,i) + ((XB(:,i) - meanB)' * V_B(:,k)) * V_B(:,k);
    end
end

disp('so far so good 4')

%% reshaping the new blocks again and putting everything back together

% create J, the approximated image file 
J = zeros(768, 1024, 3); % initialising J

for k = 1:12288
    % i and j are the starting indices of image matrix
    % where we have to fit our current block
    i = fix((k-1)/128);
    j = mod(k-1,128);
    
    % x is the 64x1 column vector which is reshaped to 8x8 block size
    x = JR(:,k);
    x = reshape(x, [8, 8]);
    
    % Now fit that block into the new image J in the red channel
    J(8*i+1:8*(i+1), 8*j+1:8*(j+1),1) = x;
    
    % Similarly for green and blue channels
    x = JG(:,k);
    x = reshape(x, [8, 8]);
    J(8*i+1:8*(i+1), 8*j+1:8*(j+1),2) = x;
    
    x = JB(:,k);
    x = reshape(x, [8, 8]);
    J(8*i+1:8*(i+1), 8*j+1:8*(j+1),3) = x;
    
end

disp('so far so good 5')
%% view the original and compressed image
figure; imshow(I); title('Original image');
figure; imshow(J); title('Compressed image, K=5');

%% now repeat the above algorithm for K = 1...64 and plot the error
% including less comments as the code was previously explained above
% This section of code can take a long time to execute on the computer
% especially for higher values of K

% array to store the errors for each K 
error = 1:64; % initialising with 64 elements

for K = 1:64
    % computing the compressed image blocks from principal eigenvectors
    JR = zeros(64,12288);
    JG = zeros(64,12288);
    JB = zeros(64,12288);
    
    for i = 1:12288
        JR(:,i) = meanR;
        JG(:,i) = meanG;
        JB(:,i) = meanB;
        
        for k = 1:K
            JR(:,i) = JR(:,i) + ((XR(:,i) - meanR)' * V_R(:,k)) * V_R(:,k);
            JG(:,i) = JG(:,i) + ((XG(:,i) - meanG)' * V_G(:,k)) * V_G(:,k);
            JB(:,i) = JB(:,i) + ((XB(:,i) - meanB)' * V_B(:,k)) * V_B(:,k);
        end
    end
    
    J = zeros(768, 1024, 3); % initialising J
    
    % fitting the image blocks to form the compressed image
    for k = 1:12288
        i = fix((k-1)/128);
        j = mod(k-1,128);
        
        % red channel
        x = JR(:,k);
        x = reshape(x, [8, 8]);
        J(8*i+1:8*(i+1), 8*j+1:8*(j+1),1) = x;
        
        % blue channel
        x = JG(:,k);
        x = reshape(x, [8, 8]);
        J(8*i+1:8*(i+1), 8*j+1:8*(j+1),2) = x;
        
        % green channel
        x = JB(:,k);
        x = reshape(x, [8, 8]);
        J(8*i+1:8*(i+1), 8*j+1:8*(j+1),3) = x;
    end   
    
    % compute and store the error (frobenius norm) between I and J 
    error(K) = norm(I(:,:,1)- J(:,:,1),'fro') + norm(I(:,:,2)- J(:,:,2),'fro') + norm(I(:,:,3)- J(:,:,3),'fro'); 
    % (computes the frobenius norm error between I and J for all 3 color channels and sums them up)
    fprintf('Computed error for K = %d', K)
end 

%% plot the errors
X = 1:64;
plot(X, error)

% Comments
% From the error plot, it is evident that the error between I and J is
% monotonically decreasing (similar to an negative exponential curve)
% which helps in compressing large images with smaller values of K