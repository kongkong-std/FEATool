load('mat.mat')    % K
load('rhs.mat')    % L

d = diag(K);
d = abs(d);
d(d == 0) = 1;

d_sqrt = sqrt(d);
D_sqrt = spdiags(d_sqrt, 0, size(K, 1), size(K, 1));

d_inv_sqrt = 1 ./ sqrt(d);
D_inv_sqrt = spdiags(d_inv_sqrt, 0, size(K, 1), size(K, 1));

K = D_inv_sqrt * K * D_inv_sqrt;
L = D_inv_sqrt * L;

save('mat_scale.mat', 'K');
save('rhs_scale.mat', 'L');
save('diag_sqrt.mat', 'D_sqrt');
