A = zeros(n-1,n);
D = zeros(n-1,n);
for i = 1:n-1
    A(i,i) = 1; A(i,i+1) = 1;
    D(i,i) = 1; D(i,i+1) = -1;
end
Ji = mi*li^2/3;
J = ones(n,1)*Ji;
J = diag(J);
M = ones(n,1)*mi;
M = diag(M);
L = ones(n,1)*li;
L = diag(L);
H = L*A'*inv(D*inv(M)*D')*A*L;
N = inv(M)*D'*inv(D*inv(M)*D')*A*L;
m = n*mi;
e = ones(n,1);
e_z = zeros(n,1);
E = [e e_z;e_z e];

Ct = ones(n,1)*cti;
Ct = diag(Ct);
Cn = ones(n,1)*cni;
Cn = diag(Cn);

Df = [Ct*M zeros(size(Ct*M));zeros(size(Ct*M)) Cn*M];
Dtau = Cn*J;