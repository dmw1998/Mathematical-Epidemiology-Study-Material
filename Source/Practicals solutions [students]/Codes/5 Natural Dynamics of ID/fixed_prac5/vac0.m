function y = vac0(t)

nt = length(t);
y = zeros(nt,1);

for k = 1:nt
    if t > 18250
        y(k) = 0;
    else 
        y(k) = 0;
    end
end
end