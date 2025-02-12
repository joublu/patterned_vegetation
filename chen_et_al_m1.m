N=11;
L=22.839;
x=linspace(0,L,N);

delta=0.05;
sigma0=0.0001;
m = 1;
eps=0.005;
a0=3;
a1=1.3;
dt=0.1;
T=(a0 -a1)/eps;
tt=0:dt:T;
aa= a0 -eps*tt;

n0=a0/2+sqrt(a0^2/4-1);
n=(x') *0+n0; w=1./n;
dx=x(2)-x(1);
tout=0;

sd_list = [];
for k_int=floor(-N/2+1):floor(N/2)
    k = k_int * (L/(2 * pi));

    n_plus = @(s) ((a0 -s) / 2) + sqrt(((a0 -s)^2) / 4 - 1);
    n_plus_prime = @(s) - (a0/2 - s/2)/(2*((a0 - s)^2/4 - 1)^(1/2)) - 1/2;
    alpha = @(s) - (k^2 * delta) + 1 - (2 * (n_plus(s)^2)) / (k^2 + 1 + n_plus(s)^2);
    beta = @(s) (sigma0 * n_plus(s)^2) / (k^2 + 1 + n_plus(s)^2);
    % alpha_prime = @(s) -(4*n_plus(s)*n_plus_prime(s) * (k^2 + 1 + n_plus(s)^2) - 2*n_plus(s)*n_plus_prime(s) * 2*n_plus(s)) / (k^2 + 1 + n_plus(s)^2)^2; 
    alpha_prime = @(s) (8*(a0/2 - s/2 + ((a0 - s)^2 - 4)^(1/2)/2)*(a0 - s + ((a0 - s)^2 - 4)^(1/2)))/(((a0 - s)^2 - 4)^(1/2)*((a0 - s + ((a0 - s)^2 - 4)^(1/2))^2 + 4*k^2 + 4)) - (4*(a0 - s + ((a0 - s)^2 - 4)^(1/2))^4)/(((a0 - s)^2 - 4)^(1/2)*((a0 - s + ((a0 - s)^2 - 4)^(1/2))^2 + 4*k^2 + 4)^2);

    sp_initial_guess = 0; % always real for initial guess 
    try
        sp = fzero(alpha, sp_initial_guess);
        if ~isreal(sp)
            % warning('sp is complex for k_int = %d, skipping this iteration.', k_int);
            continue;
        end
    catch ME
        % warning(ME.identifier, 'fzero for sp failed: %s. Skipping iteration for k_int = %d.', ME.message, k_int);
        continue;
    end

    equation = @(sd) integral(alpha, sp, sd) + eps * log(beta(sp) * (pi / (eps * alpha_prime(sp)))^(1/4));

    sd_initial_guess = 180*eps;
    try
        sd = fzero(equation, sd_initial_guess);
        if ~isreal(sd)
            % warning('sd is complex for k_int = %d, skipping this iteration.', k_int);
            continue;
        end
    catch ME
        % warning('fzero for sd failed: %s. Skipping iteration for k_int = %d.', ME.message, k_int);
        continue;
    end

    sd_list(end+1) = sd;
end

disp("les valeurs dans sd_list")
disp(sd_list);


% laplacian with neumann boundary conditions
% Lap=-2*diag(ones(1,N))+diag(ones(1,N-1) ,1)+diag(ones(1,N-1) ,-1);
% Lap(1,2)=2;Lap(N,N-1)=2;

% M1=delta*Lap/dx^2+eye(N)*(-1/dt -m);
% M2=Lap/dx^2-eye(N);

% maxn=[]; spread=[];

% for idx=1:numel(tt)
%     a=aa(idx);
%     t=tt(idx);

%     noise= randn(N,1)*sqrt(dt)*sigma0*sqrt(N);
%     wnext =(M2 -diag(n.^2))\(-a-noise/dt);
%     nnext=M1\(-n/dt -n.^2.*w);
%     n=nnext; w=wnext;

%     maxn(end+1)=max(n);
%     spread(end+1) = (max(n)-min(n))/mean(n);
%     if t>tout
%         tout=tout+10;
%         figure(1);
%         plot(x,n,x,w);
%         legend(' n' ,' w' ); xlabel( ' x' );
%         title(sprintf( ' t=%g a=%g ' ,t, a));
%         drawnow;
%     end ;
% end ;

% figure(2); hold on;
% plot(aa,maxn);
% % spread contains the sample points, aa contains the corresponding values
% % value of a when the spread is equal to 1
% at = interp1( spread , aa , 1); % revoir ca
% plot([at,at],[0,3], ' --b' );
% xlabel(' a' ); ylabel( ' max(n)' );
% title(sprintf('\\delta = %g, m = %g, \\sigma_0 = %g, \\epsilon = %g',delta, m, sigma0, eps));
% text(at, 0, 'a_d', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 12, 'Color', 'blue');