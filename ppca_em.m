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

q = 2;

% init
W = ones(d, q);
sigma = 1;
epsilon = 0.01;

% loop
while (true)
    M = W'*W + sigma^2 * eye(q);
    W_new = S*W*inv(sigma^2 * eye(q) + inv(M)*W'*S*W);
    sigma_new = sqrt(1/d * trace(S - S*W*inv(M)*W_new'));
    if(abs(sigma_new - sigma) < epsilon && max(max(abs(W_new - W))) < epsilon)
        break;
    endif
    W = W_new;
    sigma = sigma_new;
endwhile

W = W_new;
sigma = sigma_new;



M = W'*W + sigma^2 * eye(q);
for i = 1:N
    Tnorm(i,:) = T(i,:) - mu;
endfor

X = inv(M) * W' * Tnorm';

plot(X(1,:), X(2,:), "linestyle", 'none');
for i = 1:N
    text(X(1,i), X(2,i), num2str(i));
endfor
