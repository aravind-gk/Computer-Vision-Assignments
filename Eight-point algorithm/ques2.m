clc
clear all
matched_pt = importdata('matches.txt');

K1 = importdata('calibration_cam1.txt');
K2 = importdata('calibration_cam2.txt');
motion = importdata('motion.txt');
N= size(matched_pt,1);
k=10;
T=15;
sample_pole_1 = zeros(T,2);
sample_pole_2 = zeros(T,2);

sample_pole_bt1 = zeros(T,2);
sample_pole_bt2 = zeros(T,2);
for i = 1:T
    y = randi([1 N],1,N-k);
    m_pt = matched_pt( y ,:);
    
    m_pt_1 = m_pt(:,1:2);
    m_pt_2 = m_pt(:,3:4);
    
    x1 = m_pt_1(:, 1);
    y1 = m_pt_1(:, 2);
    
    x2 = m_pt_2(:, 1);
    y2 = m_pt_2(:, 2);
    
    A = [x2.* x1 , x2.* y1 , x2 , y2.* x1 , y2.* y1 , y2 , x1 , y1 , ones(N-k,1)];
    [~, ~, V] = svd(A);
    
    fMatrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];
    [U, S, V] = svd(fMatrix);
    fMatrix = U(:, 1) * S(1,1) * transpose(V(:, 1)) + U(:, 2) * S(2,2) * transpose(V(:, 2));
    
    fMatrix = fMatrix/fMatrix(3,3);
    
    [epipole1,epipole2] = find_epipole(fMatrix);
    
    [~,ep_bt1] = isEpipoleInImage(fMatrix,size(im1));
    [~,ep_bt2] = isEpipoleInImage(fMatrix',size(im1));
    
    sample_pole_1(i,:)= epipole1';
    sample_pole_2(i,:)= epipole2';
    
    sample_pole_bt1(i,:) = ep_bt1;
    sample_pole_bt2(i,:) = ep_bt2;
    
end

[mean1 , cov1] = mean_cov(sample_pole_1);
[mean2 , cov2] = mean_cov(sample_pole_2);

sample_pole_hartley_1 = zeros(T,2);
sample_pole_hartley_2 = zeros(T,2);
for i = 1:T
    y = randi([1 N],1,N-k);
    m_pt = matched_pt( y ,:);
    
    m_pt_1 = m_pt(:,1:2);
    m_pt_2 = m_pt(:,3:4);
    
    
    p1 = transpose([m_pt_1, ones(N-k, 1)]);
    p2 = transpose([m_pt_2, ones(N-k, 1)]);
    norm1 = getNormMat2d(p1);
    norm2 = getNormMat2d(p2);
    
    % Normalisation
    p1 = norm1 * p1;
    p2 = norm2 * p2;
    
    p1 = transpose(p1 ./ repmat(p1(3, :), [3, 1]));
    p2 = transpose(p2 ./ repmat(p2(3, :), [3, 1]));

    x1 = p1(:, 1);
    y1 = p1(:, 2);
    x2 = p2(:, 1);
    y2 = p2(:, 2);
    
    A = [x2.* x1 , x2.* y1 , x2 , y2.* x1 , y2.* y1 , y2 , x1 , y1 , ones(N-k,1)];
    [~, ~, V] = svd(A);
    
    fMatrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];
    [U, S, V] = svd(fMatrix);
    fMatrix = U(:, 1) * S(1,1) * transpose(V(:, 1)) + U(:, 2) * S(2,2) * transpose(V(:, 2));
    Hartly_fMatrix = norm2' * fMatrix * norm1;
    
    %Hartly_fMatrix = Hartly_fMatrix/Hartly_fMatrix(3,3);
    
    [epipole1,epipole2] = find_epipole(Hartly_fMatrix);
    
    sample_pole_hartley_1(i,:)= epipole1(1:2,:)';
    sample_pole_hartley_2(i,:)= epipole2(1:2,:)'; 
end
[mean_hartley_1 , cov_hartley_1] = mean_cov(sample_pole_hartley_1);
[mean_hartley_2 , cov_hartley_2] = mean_cov(sample_pole_hartley_2);
    
    
    
    
    
    

