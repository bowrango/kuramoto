
n = 8;

rng default

x = rand(n,1);
J = rand(n);
J = J + J.';
J(1:n+1:end) = 0;

E = lyapunov(x, 1, 1, J);

expE = 0;
for ii = 1:n
    for jj = 1:n
        if ii~=jj
            expE = expE - J(ii,jj)*cos(pi*(x(ii)-x(jj)));
        end
    end
    expE = expE - cos(2*pi*x(ii));
end

TOL = 4*n*eps;
assert(norm(E-expE) < TOL)