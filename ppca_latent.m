function X = ppca_latent(T, W, sigma)
    [N, d] = size(T);
    [~, q] = size(W);

    M = W'*W + sigma^2 * eye(q);
    for j = 1:d
        mu(j) = mean(T(:,j));
    endfor
    for i = 1:N
        Tnorm(i,:) = T(i,:) - mu;
    endfor

    X = inv(M) * W' * Tnorm';
endfunction
