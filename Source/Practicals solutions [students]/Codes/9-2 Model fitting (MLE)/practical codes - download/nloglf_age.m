function val = nloglf_age(theta,uk)

Pa = uk.Pa;
Na = uk.Na;
time_stamp = uk.time_stamp;

check_age = sum(time_stamp <= 15);
before_15 = -sum(log(binopdf(Pa(1:check_age),Na(1:check_age),...
    1-exp(-theta(1)*time_stamp(1:check_age)))));
after_15 = -sum(log(binopdf(Pa(check_age+1:end),Na(check_age+1:end),...
    1-exp(-15*(theta(1)-theta(2))).*exp(-theta(2)*time_stamp(check_age+1:end)))));
val = before_15 + after_15;
