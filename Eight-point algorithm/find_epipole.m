function [epipole1,epipole2] = find_epipole(F)
[~, ~, V1] = svd(F);
[~, ~, V2] = svd(F');
epipole1 = V1(:,3)/V1(3,3);
epipole2 = V2(:,3)/V2(3,3);
epipole1 = epipole1(1:2,:);
epipole2 = epipole2(1:2,:);
end

