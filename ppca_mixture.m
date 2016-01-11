function ppca_mixture
    filename = 'data/virus3.dat';
    T = importdata(filename);
    [N, d] = size(T);
    M = 3;  % number of ppca analysers considered
    q = 2;  % dimension of ppcas

    % init
    % initializing the posteriors (p(i|t_n), indexed by (n, i). R in text)
    for n=1:N
        k = randi(M);
        for i=1:M
            if(i==k)
                posteriors(n, i) = 1;
            else
                posteriors(n, i) = 0;
            endif
        endfor
    endfor
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
            endfor
            new_mus(i, :) /= sum(posteriors(:, i));
        endfor

        % computing covariance matrices (S in text)
        covariances = cell(M);
        for i=1:M
            covariances{i} = zeros(d, d);
            for n=1:N
                covariances{i} += posteriors(n, i)*(T(n, :)' - new_mus(i, :)')*(T(n, :)' - new_mus(i, :)')';
            endfor
            covariances{i} = covariances{i} / (new_priors(i) * N);
        endfor


        % applying ppca using covariance matrices
        % (Ws are weight matrices; sigmas are variances)
        new_Ws = cell(M);
        for i=1:M
            [new_Ws{i}, new_sigmas(i)] = ppca_from_covariance(covariances{i}, q);
        endfor

        % convergence check
        if(entered && max(abs(new_priors - priors)) < epsilon && max(max(abs(new_mus - mus))) < epsilon && max(max(max(abs(cell2mat(new_Ws) - cell2mat(Ws))))) < epsilon && max(abs(new_sigmas - sigmas)) < epsilon)
            break;
        endif

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
            endfor
        endfor
        % computing the joint probabilities (p(t_n, i))
        for i=1:M
            joint(:, i) = likelihoods(:, i) * priors(i);
        endfor
        % computing the data priors (p(t_n))
        data_priors = sum(joint');
        % computing and displaying the likelihood of the data set
        disp(sum(log(data_priors)));
        % computing the posteriors (p(i|t_n), indexed by (n, i). R in text)
        for n=1:N
            posteriors(n, :) = joint(n, :) / data_priors(n);
        endfor

        % we went through the loop at least once
        entered = true;
    endwhile


    % computing the latent variables
    latent = cell(M);
    for i=1:M
        latent{i} = ppca_latent(T, Ws{i}, sigmas(i));
    endfor
    % selecting likely points for each ppca analyser
    for n=1:N
        [~, i] = sort(posteriors(n, :), 'descend');
        points_to_classes(n) = i(1);
    endfor
    classes_to_points = cell(M);
    classes_to_point_numbers = cell(M);
    for n=1:N
        classes_to_points{points_to_classes(n)}(:, end+1) = latent{points_to_classes(n)}(:, n);
        classes_to_point_numbers{points_to_classes(n)}(end+1) = n;
    endfor

    % display the results of the automatic classification
    for i=1:M
        disp(classes_to_point_numbers{i});
    endfor
endfunction
