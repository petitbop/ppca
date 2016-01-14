function [W, sigma] = ppca_from_data(T, q)
    [N, d] = size(T);
    
    for j = 1:d
        mu(j) = mean(T(:,j));
    end
    S = zeros(d);
    for n = 1:N
        S = S + (T(n,:)' - mu') * (T(n,:)' - mu')';
    end
    S = 1/N * S;

    [W, sigma] = ppca_from_covariance(S, q);
end
