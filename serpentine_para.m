%% author Lzh
n = 10;
mi = 1;
li = 8;%cm
cti = 0.1;
cni = 10;
a = pi/2;
b = 2*pi;
c = 0;
% syms c real;
omega = 0.15/3;
syms theta [n 1] real;
syms t real;
for i = 1:n
    theta(i) = a*cos(i*b/n+omega*b*t)+i*c/n+omega*c*t;
end

serpentine_matrix

phi = simplify(D*theta);
phi_dot = diff(phi);

Sth = diag(sin(theta));
Cth = diag(cos(theta));
omegath = [Cth -Sth;Sth Cth];
J_t = J+Sth*H*Sth+Cth*H*Cth;
C_t = Sth*H*Sth-Cth*H*Cth;
L_t = [Sth*N' -Cth*N']';