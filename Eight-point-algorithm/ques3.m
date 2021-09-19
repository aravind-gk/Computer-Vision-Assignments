im1 = im2single(imread('Gandhi_1.jpg'));
im2 = im2single(imread('Gandhi_2.jpg'));
im1g = rgb2gray(im1);
im2g = rgb2gray(im2);
[f1,d1] = vl_sift(im1g) ;
[f2,d2] = vl_sift(im2g) ;
[matches, ~] = vl_ubcmatch(d1,d2);

X1 = f1(1:2,matches(1,:)) ; X1(3,:) = 1 ;
X1 = X1';
X2 = f2(1:2,matches(2,:)) ; X2(3,:) = 1 ;
X2 = X2';
N = size(matches,2) ;
K1 = importdata('calibration_cam1.txt');
K2 = importdata('calibration_cam2.txt');
motion = importdata('motion.txt');
Rg = motion(:,1:3);

Max_fMatrix = zeros(3);
Max_hartley_fMatrix = zeros(3);
Max_score=-999;
Max_hartley_score = -999;

for i= 1:1000
    
    y = randi([1 N],1,8);
    m_pt_1 = X1(y,1:2);
    m_pt_2 = X2(y,1:2);
    [E,HE,F,HF] = eightPoint(m_pt_1, m_pt_2, K1, K2);
    score = cal_score(X1,X2,F);
    score_hartley = cal_score(X1,X2,HF);
    
    if score > Max_score
        Max_score = score;
        in1=m_pt_1;
        Max_fMatrix = E;
    end
    
    if score_hartley > Max_hartley_score
        in2 = m_pt_2;
        Max_hartley_score = score_hartley ;
        Max_hartley_fMatrix = HE;
    end
end

[Re,te] = Rotation_and_trans(Max_fMatrix);
[h_Re , h_te] = Rotation_and_trans(Max_hartley_fMatrix);

[theta,n] = theta_and_axis(Re * inv(Rg));
[h_theta,h_n] = theta_and_axis(h_Re * inv(Rg));
theta
h_theta





















