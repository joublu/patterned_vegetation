syms s a0 k delta real

n_plus_sym = (a0 - s)/2 + sqrt(((a0 - s)^2)/4 - 1);
alpha_sym = - (k^2 * delta) + 1 - (2 * (n_plus_sym^2)) / (k^2 + 1 + n_plus_sym^2);

% Compute the derivative of n_plus with respect to s
n_plus_prime_sym = simplify(diff(n_plus_sym, s));
alpha_prime_sym = simplify(diff(alpha_sym, s));

% Display the result
disp('The derivative of n_plus(s) is:');
% disp(n_plus_prime_sym)
pretty(alpha_prime_sym)