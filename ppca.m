clear

filename = 'data/virus3.dat';
T = importdata(filename);
[N, d] = size(T);

for j = 1:d
    mu(j) = mean(T(:,j));
endfor

S = zeros(d);
for n = 1:N
    S = S + (T(n,:)' - mu') * (T(n,:)' - mu')';
endfor
S = 1/N * S;

[Wpca, lambda] = eig(S);
lambda = diag(lambda);
[lambda, i] = sort(lambda, 'descend');
Wpca = Wpca(:,i);
q = 2;

U = Wpca(:,1:q);
L = diag(lambda)(1:q, 1:q);
sigma = sqrt(1/(d-q) * sum(lambda(q+1:d)));
W = U * sqrt(L - sigma^2*eye(q));

M = W'*W + sigma^2 * eye(q);
for i = 1:N
    Tnorm(i,:) = T(i,:) - mu;
endfor

X = inv(M) * W' * Tnorm';
%disp(X);
plot(X(1,:), X(2,:), "linestyle", 'none');
for i = 1:N
    text(X(1,i), X(2,i), num2str(i));
endfor
