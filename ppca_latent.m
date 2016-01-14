function X = ppca_latent(T, W, sigma)
    [N, d] = size(T);
    [~, q] = size(W);

    M = W'*W + sigma^2 * eye(q);
    for j = 1:d
        mu(j) = mean(T(:,j));
    end
    for i = 1:N
        Tnorm(i,:) = T(i,:) - mu;
    end

    X = inv(M) * W' * Tnorm';
end
