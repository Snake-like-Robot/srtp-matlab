syms t real;
phi = cell(n-1,1);
for i = 1:n-1
    phi{i} = alpha*sin(omega*t+(i-1)*beta)+gamma;
end
phi_dot = cellfun(@diff,phi);