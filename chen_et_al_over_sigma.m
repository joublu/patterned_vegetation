N = 100;
L = 22.839;
x = linspace(0, L, N);
delta = 0.04;
m = 1; %new parameter m
a0 = 2.5;
a1 = 1.2;
dt = 0.1;
eps = 0.025;

sigma0_values = 0:0.0001:0.0005;  
num_sigma = numel(sigma0_values);

colors = interp1([0, 0.5, 1], [0, 1, 0; 1, 1, 0; 1, 0, 0], linspace(0, 1, num_sigma));

figure(2); hold on;

for i = 1:num_sigma
    sigma0 = sigma0_values(num_sigma-i+1);
    
    T = (a0 - a1) / eps;
    tt = 0:dt:T;
    aa = a0 - eps * tt;
    
    n0 = a0/2 + sqrt(a0^2/4 - 1);
    n = repmat(n0, [N, 1]);
    w = 1./n;
    dx = x(2)-x(1);
    tout = 0;
    
    % Laplacian with Neumann boundary conditions
    Lap = -2 * diag(ones(1, N)) + diag(ones(1, N-1),1) + diag(ones(1, N-1),-1);
    Lap(1,2) = 2;
    Lap(N, N-1) = 2;
    
    M1 = delta * Lap/dx^2 + eye(N)*(-1/dt - m);
    M2 = Lap/dx^2 - eye(N);
    
    maxn = [];
    spread = [];
    
    for idx = 1:numel(tt)
        a = aa(idx);
        t = tt(idx);
        
        noise = randn(N, 1)*sqrt(dt)*sigma0*sqrt(N);
        wnext = (M2 - diag(n.^2)) \ (-a - noise/dt);
        nnext = M1 \ (-n/dt - n.^2.*w);
        n = nnext;
        w = wnext;
        
        maxn(end+1) = max(n);
        spread(end+1) = (max(n) - min(n))/mean(n);
    end
    h = plot(aa, maxn, 'Color', colors(i, :), 'DisplayName', sprintf('\\sigma_0 = %g', sigma0));
    try
        at = interp1(spread, aa, 1);
        plot([at, at], [0, 3], '--', 'Color', h.Color, 'HandleVisibility', 'off');
    catch
        % warning('a_d not found for \\sigma_0 = %g', sigma0);
    end
end

xlabel('a');
ylabel('max(n)');
title(sprintf('\\delta = %g, m = %g, \\epsilon = %g', delta, m, eps));
legend('Location','northeastoutside');

hold on;
box on;
a = a1:0.001:a0;
Ep_n = @(a) a./(2*m) - sqrt((a./(2*m)).^2 - 1);
Em_n = @(a) a./(2*m) + sqrt((a./(2*m)).^2 - 1);
plot(a, Ep_n(a), 'b:', 'DisplayName', 'E_{+}');
plot(a, Em_n(a), 'b', 'DisplayName', 'E_{-}');
xlabel('a');
ylabel('n_{eq}');
% plot([2*m 2*m], [0 1], 'k:', 'HandleVisibility', 'off'); % a_C

legend show;
hold off;