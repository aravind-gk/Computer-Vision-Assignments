%function mosaic = merge(file1, file2)
function mosaic = merge(image_1, image_2)

% convert image to single precision [0-255 converted to 0-1]
image_1 = im2single(image_1) ;
image_2 = im2single(image_2) ;

% convert the given images to grayscale so that SIFT function can be used
% in one go, instead of calling SIFT for each color channel

image_1_grey = rgb2gray(image_1);
image_2_grey = rgb2gray(image_2);

% Use SIFT library to compute features and descriptors from both images
% f are the features, of size 4xN
% d are descriptors, of size 128xN

[f1,d1] = vl_sift(image_1_grey) ;
[f2,d2] = vl_sift(image_2_grey) ;

% find and match the best features accross both images according to
% best/minimum euclidean distance stored in 'scores' 
[matches, scores] = vl_ubcmatch(d1,d2);

% find the total number of feature matches 
% here '2' denotes columns 
numMatches = size(matches,2) ;

% f = [x,y,scale,orientation]
% select x and y from f, and only those corresponding to matches 
% set third coordinate as 1 to make homogenous coordinate
X1 = f1(1:2,matches(1,:)) ; X1(3,:) = 1 ;
X2 = f2(1:2,matches(2,:)) ; X2(3,:) = 1 ;

%% Compute the homography matrix using RANSAC algorithm

% clear the workspace variables after each function call
clear H score ok ;

for t = 1:100
  % estimate homography matrix using any 4 random points
  % select 4 random matching features among all matches
  subset = vl_colsubset(1:numMatches, 4) ;
  A = [] ;
  for i = subset
      % find H matrix using method from class
    A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
  end
  % Compute SVD of A to find the eigenvector of A'A corresponding to the
  % smallest eigenvalue of A'A, which is the solution of the linear
  % system of equations given any 4 features
  [U,S,V] = svd(A) ;
  % V(:,9) pick the last eigenvector in the 12x9 matrix V
  % { } dynamically append element to list
  H{t} = reshape(V(:,9),3,3) ;

  % score homography
  X2_ = H{t} * X1 ; % points after applying the transformation H
  % check whether the line obtained from RANSAC satisfy the threshold value
  du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
  dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
  % count the number of 1's in the binary vector ok (which indicate the
  % number of inlier points, and compute a score to determine whether this
  % is a good prediction or not
  ok{t} = (du.*du + dv.*dv) < 6*6 ;
  score(t) = sum(ok{t}) ;
end
% choose the line which has the maximum score given by RANSAC (in the form
% of max number of inlier points)
[score, best] = max(score) ;
H = H{best} ;
ok = ok{best} ;

%% Create a mosaic from the given features and homography matrix



box2 = [1  size(image_2,2) size(image_2,2)  1 ;
        1  1           size(image_2,1)  size(image_2,1) ;
        1  1           1            1 ] ;
%box2_ = inv(H) * box2 ;
box2_ = H \ box2;
% ./ is element wise division in a vector (in matlab) (dividing by the
% scaling factor)
box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;

% x_axis_endpoints stores the range of x-coordinates of box2_ (which defines its margins)
% similarly y_axis_endpoints stores the y-coordinates of box2_ from min to max
x_axis_endpoints = min([1 box2_(1,:)]):max([size(image_1,2) box2_(1,:)]) ;
y_axis_endpoints = min([1 box2_(2,:)]):max([size(image_1,1) box2_(2,:)]) ;

% create a meshgrid for box2_ 
[u,v] = meshgrid(x_axis_endpoints,y_axis_endpoints);
% perform backward warping of image1 on canvas box2_
% im2double converts pixel values from 0...1 to 0...255

im1_ = vl_imwbackward(im2double(image_1),u,v) ;

% u_, v_, z_ are the homogenous coordinates of second image
z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;

im2_ = vl_imwbackward(im2double(image_2),u_,v_) ;

% mass is 2 for those areas of the image which are overlapping (in this
% case the value of mass is 2, hence dividing the sum of pixels by mass
% negates this overlapping). In non-overlappping portions of the image, the
% mass is 1, hence dividing by the mass for non-overlapping pixels makes no
% difference
% note that for overlapping pixels, we are taking the average of the two
% values

%% No blending (Comment out this section and comment below section to get non-blending output)

%mosaic = im1_;
%mosaic(isnan(mosaic)) = im2_(isnan(mosaic));
%% Average Blending

mass = ~isnan(im1_) + ~isnan(im2_) ;
im1_(isnan(im1_)) = 0 ;
im2_(isnan(im2_)) = 0 ;
mosaic = (im1_ + im2_) ./ mass ; % blending to average out the overlapping pixels 

end