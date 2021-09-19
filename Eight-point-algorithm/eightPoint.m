function[eMatrix, Hartly_eMatrix, fMatrix , Hartly_fMatrix] = eightPoint(m_pt_1, m_pt_2, K1, K2)
[n1, c1] = size(m_pt_1);
[n2, c2] = size(m_pt_2);
if((c1 ~= 2) || (c2 ~= 2))
    error('Points are not formated with correct number of coordinates.');
end
if((n1 < 8) || (n2 < 8))
    error('There are not enough points to carry out the operation.');
end
p1 = m_pt_1(1:8,:);
p2 = m_pt_2(1:8,:);

x1 = p1(:, 1);
y1 = p1(:, 2);

x2 = p2(:, 1);
y2 = p2(:, 2);

A = [x2.* x1 , x2.* y1 , x2 , y2.* x1 , y2.* y1 , y2 , x1 , y1 , ones(8,1)];

[~, ~, V] = svd(A);
fMatrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];

[U, S, V] = svd(fMatrix);
fMatrix = U(:, 1) * S(1,1) * transpose(V(:, 1)) + U(:, 2) * S(2,2) * transpose(V(:, 2));
eMatrix = K2' * fMatrix * K1;

p1 = transpose([m_pt_1(1: 8, :), ones(8, 1)]);
p2 = transpose([m_pt_2(1: 8, :), ones(8, 1)]);
norm1 = getNormMat2d(p1);
norm2 = getNormMat2d(p2);

p1 = norm1 * p1;
p2 = norm2 * p2;

p1 = transpose(p1 ./ repmat(p1(3, :), [3, 1]));
p2 = transpose(p2 ./ repmat(p2(3, :), [3, 1]));

x1 = p1(:, 1);
y1 = p1(:, 2);
x2 = p2(:, 1);
y2 = p2(:, 2);

A = [x2.* x1 , x2.* y1 , x2 , y2.* x1 , y2.* y1 , y2 , x1 , y1 , ones(8,1)];
[~, ~, V] = svd(A);
Hartly_fMatrix = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];
[U, S, V] = svd(Hartly_fMatrix);

Hartly_fMatrix = U(:, 1) * S(1,1) * transpose(V(:, 1)) + U(:, 2) * S(2,2) * transpose(V(:, 2));
Hartly_fMatrix = norm2' * Hartly_fMatrix * norm1;
Hartly_eMatrix = K2' * Hartly_fMatrix * K1;

end
