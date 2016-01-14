function [W, sigma] = ppca_from_covariance(S, q)
    [d, ~] = size(S);

    [Wpca, lambda] = eig(S);
    lambda = diag(lambda);
    [lambda, i] = sort(lambda, 'descend');
    Wpca = Wpca(:,i);
    U = Wpca(:,1:q);
    L = diag(lambda)(1:q, 1:q);

    sigma = sqrt(1/(d-q) * sum(lambda(q+1:d)));
    W = U * sqrt(L - sigma^2*eye(q));
end
