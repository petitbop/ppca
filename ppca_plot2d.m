function ppca_plot2d(X)
    [~, N] = size(X);

    axis([-3, 3, -3, 3])
    for i = 1:N
        text(X(1,i), X(2,i), num2str(i));
    end
end
