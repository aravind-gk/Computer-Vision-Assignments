function Nmatrix = getNormMat2d(x)

centroid = mean(x, 2);
dist = sqrt(sum((x - repmat(centroid, 1, size(x, 2))) .^ 2, 1));
mean_dist = mean(dist);
Nmatrix = [sqrt(2) / mean_dist, 0, -sqrt(2) / mean_dist * centroid(1);...
           0, sqrt(2) / mean_dist, -sqrt(2) / mean_dist * centroid(2);...
           0, 0, 1];

end
