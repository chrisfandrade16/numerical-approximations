mu = 0.012277471;
mu_hat = 1 - mu;
first = 0;
last = 17.1;

% w = u1_prime
% w_prime = u1_prime_prime
% v = u2_prime
% v_prime = u2_prime_prime

w = @(u1, u1_prime, u2, u2_prime) u1_prime;
w_prime = @(u1, u1_prime, u2, u2_prime) u1 + (2 * u2_prime) - ((mu_hat) * ((u1 + mu) / ((((u1 + mu)^2) + u2^2)^(3/2)))) - (mu * ((u1 - mu_hat) / ((((u1 - mu_hat)^2) + u2^2)^(3/2))));
v = @(u1, u1_prime, u2, u2_prime) u2_prime;
v_prime = @(u1, u1_prime, u2, u2_prime) u2 - (2 * u1_prime) - ((mu_hat) * (u2 / ((((u1 + mu)^2) + u2^2)^(3/2)))) - (mu * (u2 / ((((u1 - mu_hat)^2) + u2^2)^(3/2))));

steps_list = {100, 1000, 10000, 20000};
for i = 1:length(steps_list)
    steps = steps_list{i};
    h = (last - first) / steps;

    u1 = zeros(1, steps);
    u1(1, 1) = 0.994;

    u1_prime = zeros(1, steps);
    u1_prime(1, 1) = 0;

    u2 = zeros(1, steps);
    u2(1, 1) = 0;

    u2_prime = zeros(1, steps);
    u2_prime(1, 1) = -2.001585106379082522420537862224;

    u1_prime_prime = zeros(1, steps);
    u1_prime_prime(1, 1) = w_prime(u1(1, 1), u1_prime(1, 1), u2(1, 1), u2_prime(1, 1));

    u2_prime_prime = zeros(1, steps);
    u2_prime_prime(1, 1) = v_prime(u1(1, 1), u1_prime(1, 1), u2(1, 1), u2_prime(1, 1));

    x = zeros(1, steps);
    x(1, 1) = first;

    n = 1;
    
    while n < steps
        % y_n
        u1_n = u1(1, n);
        u1_n_prime = u1_prime(1, n);
        u2_n = u2(1, n);
        u2_n_prime = u2_prime(1, n);

        % k1 = f(y_n)
        k1_w = h * w(u1_n, u1_n_prime, u2_n, u2_n_prime);
        k1_w_prime = h * w_prime(u1_n, u1_n_prime, u2_n, u2_n_prime);
        k1_v = h * v(u1_n, u1_n_prime, u2_n, u2_n_prime);
        k1_v_prime = h * v_prime(u1_n, u1_n_prime, u2_n, u2_n_prime);

        % y_n + h * k1/2
        u1_n_k1 = u1_n + ((k1_w / 2));
        u1_n_k1_prime = u1_n_prime + ((k1_w_prime / 2));
        u2_n_k1 = u2_n + ((k1_v / 2));
        u2_n_k1_prime = u2_n_prime + ((k1_v_prime / 2));

        % k2 = f(y_n + h * k1/2)
        k2_w = h * w(u1_n_k1, u1_n_k1_prime, u2_n_k1, u2_n_k1_prime);
        k2_w_prime = h * w_prime(u1_n_k1, u1_n_k1_prime, u2_n_k1, u2_n_k1_prime);
        k2_v = h * v(u1_n_k1, u1_n_k1_prime, u2_n_k1, u2_n_k1_prime);
        k2_v_prime = h * v_prime(u1_n_k1, u1_n_k1_prime, u2_n_k1, u2_n_k1_prime);

        % y_n + h * k2/2
        u1_n_k2 = u1_n + ((k2_w / 2));
        u1_n_k2_prime = u1_n_prime + ((k2_w_prime / 2));
        u2_n_k2 = u2_n + ((k2_v / 2));
        u2_n_k2_prime = u2_n_prime + ((k2_v_prime / 2));

        % k3 = f(y_n + h * k2/2)
        k3_w = h * w(u1_n_k2, u1_n_k2_prime, u2_n_k2, u2_n_k2_prime);
        k3_w_prime = h * w_prime(u1_n_k2, u1_n_k2_prime, u2_n_k2, u2_n_k2_prime);
        k3_v = h * v(u1_n_k2, u1_n_k2_prime, u2_n_k2, u2_n_k2_prime);
        k3_v_prime = h * v_prime(u1_n_k2, u1_n_k2_prime, u2_n_k2, u2_n_k2_prime);

        % y_n + h * k3
        u1_n_k3 = u1_n + (k3_w);
        u1_n_k3_prime = u1_n_prime + (k3_w_prime);
        u2_n_k3 = u2_n + (k3_v);
        u2_n_k3_prime = u2_n_prime + (k3_v_prime);

        % k4 = f(y_n + h * k3)
        k4_w = h * w(u1_n_k3, u1_n_k3_prime, u2_n_k3, u2_n_k3_prime);
        k4_w_prime = h * w_prime(u1_n_k3, u1_n_k3_prime, u2_n_k3, u2_n_k3_prime);
        k4_v = h * v(u1_n_k3, u1_n_k3_prime, u2_n_k3, u2_n_k3_prime);
        k4_v_prime = h * v_prime(u1_n_k3, u1_n_k3_prime, u2_n_k3, u2_n_k3_prime);

        % y_n+1 = y_n + 1/6 * (k1 + 2*k2 + 2*k3 + k4) * h
        u1(1, n + 1) = u1_n + ((1/6) * (k1_w + (2 * k2_w) + (2 * k3_w) + k4_w));
        u1_prime(1, n + 1) = u1_n_prime + ((1/6) * (k1_w_prime + (2 * k2_w_prime) + (2 * k3_w_prime) + k4_w_prime));
        u2(1, n + 1) = u2_n + ((1/6) * (k1_v + (2 * k2_v) + (2 * k3_v) + k4_v));
        u2_prime(1, n + 1) = u2_n_prime + ((1/6) * (k1_v_prime + (2 * k2_v_prime) + (2 * k3_v_prime) + k4_v_prime));

        % update u1'' and u2'' with calculated u1, u1', u2, u2'
        u1_prime_prime(1, n + 1) = w_prime(u1(1, n + 1), u1_prime(1, n + 1), u2(1, n + 1), u2_prime(1, n + 1));
        u2_prime_prime(1, n + 1) = v_prime(u1(1, n + 1), u1_prime(1, n + 1), u2(1, n + 1), u2_prime(1, n + 1));

        % t_n+1 + t_n + h
        x(1, n + 1) = x(1, n) + h;

        % next step
        n = n + 1;
    end

    figure("Name", num2str(steps, "u1 vs u2 for %d steps"));
    plot(u1, u2, "r");
    hold on;

    figure("Name", num2str(steps, "u1' vs u2' for %d steps"));
    plot(u1_prime, u2_prime, "r");
    hold on;

    figure("Name", num2str(steps, "u1'' vs u2'' for %d steps"));
    plot(u1_prime_prime, u2_prime_prime, "r");
    hold on;
end