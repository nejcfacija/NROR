function mcc_pi()
    n = 1000000;
    r = 1.0;
    tocke_v_krogu = 0;
    priblizek_pi = zeros(n, 1);

    for i = 1:n
        x = -r + 2 * r * rand();
        y = -r + 2 * r * rand();
        if x^2 + y^2 <= r^2
            tocke_v_krogu = tocke_v_krogu + 1;
        end
        priblizek_pi(i) = (tocke_v_krogu / i) * 4;
    end

    figure;
    plot(1:n, priblizek_pi);

    fprintf('Končni približek π: %f\n', priblizek_pi(n));
end

