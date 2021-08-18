function [R,t] = Rotation_and_trans(E)

W1 = [0 -1 0;1 0 0;0 0 1];
W2 = W1';
[U,~,V] = svd(E);
R = (-1)*U * W1 * V';
t = U(:,3);
t = t/norm(t);
end
