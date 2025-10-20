function E = lyapunov(x, K, Ks, J)
% Lyapunov (4.7)
shift = x - x.';
E = -K*cos(pi*(shift(:)))'*J(:) - Ks*sum(cos(2*pi*x));
end