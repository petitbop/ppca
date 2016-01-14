function [W, sigma] = em_ppca_from_covariance(S, q)
    [d, ~] = size(S);
    
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
        end
        W = W_new;
        sigma = sigma_new;
    end
    
    W = W_new;
    sigma = sigma_new;
end
