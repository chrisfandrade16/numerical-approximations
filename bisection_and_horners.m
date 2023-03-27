tol = 10^(-6);
A = 1.94;
B = 2.07;
X = linspace(A, B);
f1 = @(x) (x-2)^9;
f2 = [1, -18, 144, -672, 2016, -4032, 5376, -4608, 2304, -512];
y1 = zeros(length(X), 1);
y2 = zeros(length(X), 1);

for index=1:length(X)
   y1(index) = f1(X(index));
   y2(index) = horner_evaluate(f2, X(index));
end

plot(X, y1, "r");
hold on;
plot(X, y2, "b");
legend("f1(x)=(x-2)^9", "f2(x)=x^9-18x^8+144x^7-672x^6+2016x^5-4032x^4+5376x^3-4608x^2+2304x-512");

y1root = bisection_method(f1, A, B, tol, false);
y2root = bisection_method(f2, A, B, tol, true);

fprintf("f1 = 0 when x = %.32f, using bisection method\n", y1root);
fprintf("f2 = 0 when x = %.32f, using bisection method\n", y2root);

function result = horner_evaluate(constants, x)
    last_index = length(constants);

    if last_index == 1
        result = constants(1); 
    else
        result = constants(last_index) + (x * horner_evaluate(constants(1:last_index - 1), x));
    end
end

function root = bisection_method(f, a, b, tol, useHorner)
    c = (a + b) / 2;

    if useHorner
        f_at_a = horner_evaluate(f, a);
        f_at_c = horner_evaluate(f, c);
    else
        f_at_a = f(a);
        f_at_c = f(c);
    end

    if f_at_c == 0
        root = c;
    elseif (b - a) < tol
        root = c;
    else
        if sign(f_at_c) == sign(f_at_a)
           root = bisection_method(f, c, b, tol, useHorner);
        else
           root = bisection_method(f, a, c, tol, useHorner);
        end
    end
end

