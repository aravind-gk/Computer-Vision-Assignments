function IMGSEG(img)

    img = rgb2gray(img);        %grayscale conversion

    img = imresize(img,[100 100]);  %resizing
    [m, n]=size(img);

    img = double(img) / 255;    %normalizing to (0,1)


    % calculating similarity matrix W
    W = zeros(m*n,m*n);   
    r = 7;
    sigma_i = 0.1;
    sigma_x = 6;

    for i=1:m
        for j=1:n
            for k=1:m
                for l=1:n
                    x = sqrt((i-k)^2 + (j-l)^2);
                    if x < r
                        W((i-1)*n+j, (k-1)*n+l) = exp(-(img(i,j)-img(k,l))^2/sigma_i^2) * exp(-x^2/sigma_x^2);
                    end
                    if i*(n-1)+j == k*(n-1)+l
                        W((i-1)*n+j, (k-1)*n+l) = 1;
                    end
                end
            end
        end
    end

    % calculating diagonal matrix D
    D = zeros(m*n,m*n);
    for i=1:m*n
        D(i,i) = sum(W(i,:));
    end


    [vec,val] = eigs(D-W,D,4,'sm');     %finding 4 smallest eigen vectors

    % 1st segment
    vec1 = abs(vec(:,2));                               
    vec1 = (vec1-min(vec1)) / (max(vec1)-min(vec1));    
    output1 = zeros(m,n);
    k = 1;
    for i=1:m
        for j=1:n
            if vec1(k) > 0.5
                output1(i,j) = 1;
            else
                output1(i,j) = 0;
            end
            k = k+1;
        end
    end

    % 2nd segment
    vec2 = abs(vec(:,3));                              
    vec2 = (vec2-min(vec2)) / (max(vec2)-min(vec2));    
    output2 = zeros(m,n);
    k = 1;
    for i=1:m
        for j=1:n
            if vec2(k) > 0.5
                output2(i,j) = 1;
            else
                output2(i,j) = 0;
            end
            k = k+1;
        end
    end


    % 3rd segment
    vec3 = abs(vec(:,4));                               
    vec3 = (vec3-min(vec3)) / (max(vec3)-min(vec3));  
    output3 = zeros(m,n);
    k = 1;
    for i=1:m
        for j=1:n
            if vec3(k) > 0.5
                output3(i,j) = 1;
            else
                output3(i,j) = 0;
            end
            k = k+1;
        end
    end

%     figure;
%     subplot(1,5,1); imshow(img); title('input image');              %printing input image
%     subplot(1,5,2); imshow(output1); title('output : segment 1');   %printing 1st segment
%     subplot(1,5,3); imshow(output2); title('output : segment 2');   %printing 2nd segment
%     subplot(1,5,4); imshow(output3); title('output : segment 3');   %printing 3rd segment


    final_output = zeros(m,n,3);
    final_output(:,:,1) = output1;
    final_output(:,:,2) = output2;
    final_output(:,:,3) = output3;
%     subplot(1,5,5); imshow(final_output); title('final output');

    figure;
    subplot(1,2,1); imshow(img); title('input image');
    subplot(1,2,2); imshow(final_output); title('output segments');

end

