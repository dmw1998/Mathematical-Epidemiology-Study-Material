function y = solve_catalytic_age(theta,a)

check_age = sum(a <= 15);
y = zeros(size(a));
y(1:check_age) = 1 - exp(-theta(1)*a(1:check_age));
y(check_age+1:end) = 1 - exp(-15*(theta(1)-theta(2))).*exp(-theta(2)*a(check_age+1:end));