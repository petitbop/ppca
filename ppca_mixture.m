function ppca_mixture
    filename = 'data/virus3.dat';
    T = importdata(filename);
    [N, d] = size(T);
    M = 3;  % number of ppca analysers considered
    q = 2;  % dimension of ppcas

    % init
    % initializing the posteriors (p(i|t_n), indexed by (n, i). R in text)
    classes = [3, 3, 3, 3, 2, 2, 3, 1, 3, 3, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2];
    for n=1:N
        for i=1:M
            if(i==classes(n))
                posteriors(n, i) = 1;
            else
                posteriors(n, i) = 0;
            end
        end
    end
    % precision for convergence checking
    epsilon = 0.01;
    entered = false;

    % loop
    while(true)
        % updating priors (p(i), pi in text) and mean vectors (mu in text)
        new_priors = 1/N * sum(posteriors);
        new_mus = zeros(M, d);
        for i=1:M
            for n=1:N
                new_mus(i, :) += posteriors(n, i) * T(n, :);
            end
            new_mus(i, :) /= sum(posteriors(:, i));
        end

        % computing covariance matrices (S in text)
        covariances = cell(M);
        for i=1:M
            covariances{i} = zeros(d, d);
            for n=1:N
                covariances{i} += posteriors(n, i)*(T(n, :)' - new_mus(i, :)')*(T(n, :)' - new_mus(i, :)')';
            end
            covariances{i} = covariances{i} / (new_priors(i) * N);
        end


        % applying ppca using covariance matrices
        % (Ws are weight matrices; sigmas are variances)
        new_Ws = cell(M);
        for i=1:M
            [new_Ws{i}, new_sigmas(i)] = ppca_from_covariance(covariances{i}, q);
        end

        % convergence check
        if(entered && max(abs(new_priors - priors)) < epsilon && max(max(abs(new_mus - mus))) < epsilon && max(max(max(abs(cell2mat(new_Ws) - cell2mat(Ws))))) < epsilon && max(abs(new_sigmas - sigmas)) < epsilon)
            break;
        end

        % replacing old parameter values
        priors = new_priors;
        mus = new_mus;
        Ws = new_Ws;
        sigmas = new_sigmas;

        % computing the likelihoods (p(t_n|i))
        for i=1:M
            C = sigmas(i)^2 * eye(d) + Ws{i}*Ws{i}';
            detC = det(C);
            invC = inv(C);
            for n=1:N
                likelihoods(n, i) = (2*pi)^(-d/2)*detC^(-1/2)*exp(-1/2*(T(n, :)' - mus(i, :)')'*invC*(T(n, :)' - mus(i, :)'));
            end
        end
        % computing the joint probabilities (p(t_n, i))
        for i=1:M
            joint(:, i) = likelihoods(:, i) * priors(i);
        end
        % computing the data priors (p(t_n))
        data_priors = sum(joint');
        % computing and displaying the likelihood of the data set
        disp(sum(log(data_priors)));
        % computing the posteriors (p(i|t_n), indexed by (n, i). R in text)
        for n=1:N
            posteriors(n, :) = joint(n, :) / data_priors(n);
        end

        % we went through the loop at least once
        entered = true;
    end


    % computing the latent variables
    latent = cell(M);
    for i=1:M
        latent{i} = ppca_latent(T, Ws{i}, sigmas(i));
    end
    % selecting likely points for each ppca analyser
    for n=1:N
        [~, i] = sort(posteriors(n, :), 'descend');
        points_to_classes(n) = i(1);
    end
    classes_to_points = cell(M);
    classes_to_point_numbers = cell(M);
    for n=1:N
        classes_to_points{points_to_classes(n)}(:, end+1) = latent{points_to_classes(n)}(:, n);
        classes_to_point_numbers{points_to_classes(n)}(end+1) = n;
    end

    % display the results of the automatic classification
    for i=1:M
        disp(classes_to_point_numbers{i});
    end

    ppca_plot2d(classes_to_points{1});
end
