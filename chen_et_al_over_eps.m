N = 20;
L = 22.839;
x = linspace(0, L, N);
delta = 0.04;
sigma0 = 0.0001;
m = 0.45;
a0=1;
a1=0.75;
dt = 0.1;
eps_values = 0.0005:0.0002:0.002;

% colormap from green to red
num_eps = numel(eps_values);
colors = interp1([0, 0.5, 1], [0, 1, 0; 1, 1, 0; 1, 0, 0], linspace(0, 1, num_eps));

figure(2); hold on;

for i = 1:num_eps
    eps = eps_values(i);
    T = (a0 - a1) / eps;
    tt = 0:dt:T;
    aa = a0 - eps * tt;
    n0=a0/2+sqrt((a0/(2*m))^2-1);
    n = (x') * 0 + n0;
    w = 1 ./ n;
    dx = x(2) - x(1);
    tout = 0;

    % Laplacian with Neumann boundary conditions
    Lap = -2 * diag(ones(1, N)) + diag(ones(1, N - 1), 1) + diag(ones(1, N - 1), -1);
    Lap(1, 2) = 2;
    Lap(N, N - 1) = 2;

    M1 = delta * Lap / dx^2 + eye(N) * (-1 / dt - m);
    M2 = Lap / dx^2 - eye(N);

    maxn = [];
    spread = [];

    for idx = 1:numel(tt)
        a = aa(idx);
        t = tt(idx);

        noise = randn(N, 1) * sqrt(dt) * sigma0 * sqrt(N);
        wnext = (M2 - diag(n.^2)) \ (-a - noise / dt);
        nnext = M1 \ (-n / dt - n.^2 .* w);
        n = nnext;
        w = wnext;

        maxn(end + 1) = max(n);
        spread(end + 1) = (max(n) - min(n)) / mean(n);
    end
    h = plot(aa, maxn, 'Color', colors(i, :), 'DisplayName', sprintf('\\epsilon = %g', eps));
    [spreadUnique, uniqueIdx] = unique(spread);
    aaUnique = aa(uniqueIdx);
    at = interp1(spreadUnique, aaUnique, 1);
    plot([at, at], [0, 3], '--', 'Color', h.Color, 'DisplayName', sprintf('a_d for \\epsilon = %g', eps), 'HandleVisibility', 'off');
end

xlabel('a');
ylabel('max(n)');
title(sprintf('\\delta = %g, m = %g, \\sigma_0 = %g', delta, m, sigma0));
legend('Location', 'northeastoutside');

% Bifurcation diagram from soresina
hold on;
box on;
a = a1:0.001:a0;
Ep_n = @(a) a./(2*m) - sqrt((a./(2*m)).^2 - 1);
Em_n = @(a) a./(2*m) + sqrt((a./(2*m)).^2 - 1);
plot(a, Ep_n(a), 'b:', 'DisplayName', 'E_{+}');
plot(a, Em_n(a), 'b', 'DisplayName', 'E_{-}');
xlabel('a');
ylabel('n_{eq}');
% axis([0.8, 1.2, 0, 2.5]);
plot([2*m 2*m], [0 1], 'k:', 'DisplayName', 'a_C'); % aC

legend show;
hold off;