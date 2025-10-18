
rng default
n = 2;
x = pi*rand(n,1);
J = rand(n);
J = J + J.';
J(1:n+1:end) = 0;

shift = x - x.';
E = sum(J(:).*cos(shift(:)));

expE = 0;
for ii = 1:n
    for jj = 1:n
        if ii~=jj
            expE = expE + J(ii,jj)*cos(x(ii)-x(jj));
        end
    end
end
norm(E-expE)