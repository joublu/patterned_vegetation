N = 20;
L = 22.839;
x = linspace(0, L, N);
y = linspace(0, L, N);
[X, Y] = meshgrid(x, y);
delta = 0.05;
sigma0 = 0.0001;
m=1;
eps = 0.005;
a0 = 2.3;
a1 = 1.3;
dt = 0.1;
T = (a0 - a1) / eps;
tt = 0:dt:T;
aa = a0 - eps * tt;
n0=a0/2+sqrt((a0/(2*m))^2-1);
n = n0 * ones(N, N);
w = 1 ./ n;
dx = x(2) - x(1);
dy = y(2) - y(1);
tout = 0;

Lap = zeros(N^2, N^2);
for i = 1:N
    for j = 1:N
        idx = (i-1)*N + j;
        Lap(idx, idx) = -4;
        i_up = (i == 1)  * N + (i ~= 1) * (i - 1); % ~= is not equal to (returns 0 or 1)
        i_down = (i == N) * 1 + (i ~= N) * (i + 1);
        j_left = (j == 1)  * N + (j ~= 1) * (j - 1);
        j_right = (j == N) * 1 + (j ~= N) * (j + 1);
        
        Lap(idx, (i_up-1)*N + j) = 1;
        Lap(idx, (i_down-1)*N + j) = 1;
        Lap(idx, (i-1)*N + j_left) = 1;
        Lap(idx, (i-1)*N + j_right) = 1;
    end
end

% disp(Lap);

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
        figure(1);
        [Xq, Yq] = meshgrid(linspace(0, L, 10*N), linspace(0, L, 4*N));
        nq = interp2(X, Y, n, Xq, Yq, 'cubic'); % to improve visuals but not mathematically correct
        imagesc(x, y, nq);
        hold on;
        % imagesc(x, y, w);
        % hold off;
        % legend('n', 'w');
        xlabel('x');
        ylabel('y');
        title(sprintf('a=%g', a));
        colorbar;
        clim([0, 3]);
        drawnow;
        % pause;
    end
end

% to plot the bifurcation diagram
figure(2); hold on;
plot(aa, maxn);
at = interp1(spread, aa, 1);
plot([at, at], [0, 3], '--b');
xlabel('a');
ylabel('max(n)');
title(sprintf('\\delta = %g, m = %g, \\sigma_0 = %g, \\epsilon = %g',delta, m, sigma0, eps));
text(at, 0, 'a_d', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'blue');