function y = vac(cov,t)

nt = length(t);
y = zeros(nt,1);

switch cov
    case 0
        y = vac0(t);
    case 60
        y = vac60(t);
    case 75
        y = vac75(t);
    case 90
        y = vac90(t);
    case 95
        y = vac95(t);
end

end