%% Given dataset (Run one section at a time)

im1 = imread('image01.jpg') ;
im2 = imread('image02.jpg') ;
im3 = imread('image03.jpg') ;
im4 = imread('image04.jpg') ;
im5 = imread('image05.jpg') ;

% Use reduced scale (0.1 or 0.2 instead of 1.0 to get faster outputs)
reduced_scale = 0.5;

im1 = imresize(im1, reduced_scale);
im2 = imresize(im2, reduced_scale);
im3 = imresize(im3, reduced_scale);
im4 = imresize(im4, reduced_scale);
im5 = imresize(im5, reduced_scale);

mosaic3 = im3;
mosaic2 = merge(mosaic3, im2);
mosaic4 = merge(mosaic2, im4);
mosaic1 = merge(mosaic4, im1);
mosaic5 = merge(mosaic1, im5);

figure(1); clf;
imagesc(mosaic5) ; 
axis image off ;
title('Resulting mosaic') ;

%% Custom dataset 
% Note that running this section different times may give different-looking
% outputs. Hence, it's best to run this code section multiple times to get the best
% output for the custom dataset.

im1 = imread('im01.jpeg') ;
im2 = imread('im02.jpeg') ;
im3 = imread('im03.jpeg') ;
im4 = imread('im04.jpeg') ;
im5 = imread('im05.jpeg') ;

% Reduce the image scale to get faster outputs
reduced_scale = 0.5;

im1 = imresize(im1, reduced_scale);
im2 = imresize(im2, reduced_scale);
im3 = imresize(im3, reduced_scale);
im4 = imresize(im4, reduced_scale);
im5 = imresize(im5, reduced_scale);

mosaic3 = im3;
mosaic2 = merge(mosaic3, im2);
mosaic4 = merge(mosaic2, im4);
mosaic1 = merge(mosaic4, im1);
mosaic5 = merge(mosaic1, im5);

figure(1); clf;
imagesc(mosaic5) ; 
axis image off ;
title('Resulting mosaic') ;
