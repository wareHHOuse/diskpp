gK = load('stiff_matrix.txt');
gM = load('mass_matrix.txt');

K = sparse(gK(:,1)+1, gK(:,2)+1, gK(:,3));
[M,N] = size(K);
M = sparse(gM(:,1)+1, gM(:,2)+1, gM(:,3), M, N);

eigs(K, M, 6, 'sm')