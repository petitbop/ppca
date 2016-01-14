function ppca_example
    filename = 'data/virus3.dat';
    T = importdata(filename);
    q = 2;

    [W, sigma] = ppca_from_data(T, q);
    X = ppca_latent(T, W, sigma);
    ppca_plot2d(X);
end
