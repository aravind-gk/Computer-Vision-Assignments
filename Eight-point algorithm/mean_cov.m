function [mean , c] = mean_cov(S)

lenth = size(S,1);

mean_x = sum(S(:,1))/lenth;
mean_y = sum(S(:,2))/lenth;
mean = [mean_x mean_y];

c = cov(S);
end

