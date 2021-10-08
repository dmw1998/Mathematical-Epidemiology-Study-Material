function [G, SG, R] = Gen1(Ini_inf, NGM, c)

% G -> NGM*Ini_inf (3,10)
% SG -> sum(G)
% R -> G_k/G_k-1   9

G = zeros(3,10);

G(:,1) = (1-c)*NGM*Ini_inf;

for i = 1:9
    G(:,i+1) = NGM*G(:,i)*(1-c);
end

SG = sum(G);
R = SG(2:end)./SG(1:end-1);

end