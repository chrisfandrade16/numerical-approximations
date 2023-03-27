function errors(tol)
    f = @(x) exp(2 .* x) .* sin(x);
    a = 0;
    b = pi/2;

    fprintf("tol=%e\n", tol);

    [composite_trapezoid_error, composite_trapezoid_n] = composite_trapezoid(f, a, b, tol);
    fprintf("trapezoid n= %d, error=%e\n", composite_trapezoid_n, composite_trapezoid_error);

    [composite_midpoint_error, composite_midpoint_n] = composite_midpoint(f, a, b, tol);
    fprintf("midpoint  n= %d, error=%e\n", composite_midpoint_n, composite_midpoint_error);

    [composite_simpson_error, composite_simpson_n] = composite_simpson(f, a, b, tol);
    fprintf("simpson   n= %d, error=%e\n", composite_simpson_n, composite_simpson_error);
end

function [composite_trapezoid_error, composite_trapezoid_n] = composite_trapezoid(f, a, b, tol)
    calculated_integral = integral(f, a, b);

    composite_trapezoid_error = 10000; % to let the while loop initially run
    composite_trapezoid_n = 0;

    r = 1;
    while composite_trapezoid_error > tol
        h = (b - a) / r;
        t = @(i) a + (i * h);

        composite_trapezoid_sum = 0;
        for i = 1:(r-1)
            composite_trapezoid_sum = composite_trapezoid_sum + f(t(i));
        end
        calculated_composite_trapezoid = ((h/2) * (f(a) + f(b))) + (h * composite_trapezoid_sum);

        composite_trapezoid_error = abs(calculated_integral - calculated_composite_trapezoid);
        composite_trapezoid_n = r;
        r = r + 1;
    end
end

function [composite_midpoint_error, composite_midpoint_n] = composite_midpoint(f, a, b, tol)
    calculated_integral = integral(f, a, b);

    composite_midpoint_error = 10000; % to let the while loop initially run
    composite_midpoint_n = 0;

    r = 1;

    while composite_midpoint_error > tol
        h = (b - a) / r;
        t = @(i) a + (i * h);

        composite_midpoint_sum = 0;
        for i = 1:r
            composite_midpoint_sum = composite_midpoint_sum + f(t(i - (1/2)));
        end
        calculated_composite_midpoint = h * composite_midpoint_sum;

        composite_midpoint_error = abs(calculated_integral - calculated_composite_midpoint);
        composite_midpoint_n = r;
        r = r + 1;
    end
end

function [composite_simpson_error, composite_simpson_n] = composite_simpson(f, a, b, tol)
    calculated_integral = integral(f, a, b);

    composite_simpson_error = 10000; % to let the while loop initially run
    composite_simpson_n = 0;

    r = 1;

    while composite_simpson_error > tol
        h = (b - a) / r;
        t = @(i) a + (i * h);

        sum_simpson_1 = 0;
        for i = 1:((r/2) - 1)
            sum_simpson_1 = sum_simpson_1 + f(t(2 * i));
        end
        sum_simpson_2 = 0;
        for i = 1:(r/2)
            sum_simpson_2 = sum_simpson_2 + f(t((2 * i) - 1));
        end
        calculated_composite_simpson = (h/3) * (f(a) + (2 * sum_simpson_1) + (4 * sum_simpson_2) + f(b));
    
        composite_simpson_error = abs(calculated_integral - calculated_composite_simpson);
        composite_simpson_n = r;
        r = r + 1;
    end
end

