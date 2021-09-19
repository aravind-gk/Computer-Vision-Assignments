function final_score = cal_score(X1,X2,fMatrix)
f = [fMatrix(1,1),fMatrix(1,2),fMatrix(1,3),fMatrix(2,1),fMatrix(2,2),fMatrix(2,3),fMatrix(3,1),fMatrix(3,2), fMatrix(3,3)]';
x1 = X1(:, 1);
y1 = X1(:, 2);

x2 = X2(:, 1);
y2 = X2(:, 2);

score = zeros(size(X1,1),1);
epi = 0.1;

A = [x2.* x1 , x2.* y1 , x2 , y2.* x1 , y2.* y1 , y2 , x1 , y1 , ones( size(X1,1),1)];
error = A*f;
score(abs(error)<epi) = 1;
final_score = sum(score);
end

