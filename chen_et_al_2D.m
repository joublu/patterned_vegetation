N = 20;
L = 22.839;
x = linspace(0, L, N);
y = linspace(0, L, N);
[X, Y] = meshgrid(x, y);
delta = 0.05;
sigma0 = 0.0001;
eps = 0.005;
a0 = 3;
a1 = 1.3;
dt = 0.1;
T = (a0 - a1) / eps;
tt = 0:dt:T;
aa = a0 - eps * tt;
n0 = a0 / 2 + sqrt(a0^2 / 4 - 1);
n = n0 * ones(N, N);
w = 1 ./ n;
dx = x(2) - x(1);
dy = y(2) - y(1);
tout = 0;

LapX = -2 * diag(ones(1, N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1);
LapX(1, 2) = 2; LapX(N, N-1) = 2;
LapY = LapX;
Lap = kron(eye(N), LapX) / dx^2 + kron(LapY, eye(N)) / dy^2;

M1 = delta * Lap + eye(N^2) * (-1 / dt - 1);
M2 = Lap - eye(N^2);

maxn = [];
spread = [];

for idx = 1:numel(tt)
    a = aa(idx);
    t = tt(idx);
    noise = randn(N, N) * sqrt(dt) * sigma0 * sqrt(N);
    noise = noise(:);
    wnext = (M2 - diag(n(:).^2)) \ (-a - noise / dt);
    nnext = M1 \ (-n(:) / dt - n(:).^2 .* w(:));
    n = reshape(nnext, N, N);
    w = reshape(wnext, N, N);
    maxn(end+1) = max(n(:));
    spread(end+1) = (max(n(:)) - min(n(:))) / mean(n(:));
    
    if t > tout
        tout = tout + 10;
        subplot(2, 1, 1);
        imagesc(x, y, n);
        hold on;
        imagesc(x, y, w);
        hold off;
        legend('n', 'w');
        xlabel('x');
        ylabel('y');
        title(sprintf('t=%g a=%g', t, a));
        colorbar;
        drawnow;
    end
end

subplot(2, 1, 2); hold on;
plot(aa, maxn);
at = interp1(spread, aa, 1);
plot([at, at], [0, 3], '--b');
xlabel('a');
ylabel('max(n)');
title(sprintf('a_d=%g', at));