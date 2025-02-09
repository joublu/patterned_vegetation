N=200;
L=22.839;
x=linspace(0,L,N);

delta=0.05;
sigma0=0; %0.0001;
m = 1; % New parameter m
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

% laplacian with neumann boundary conditions
Lap=-2*diag(ones(1,N))+diag(ones(1,N-1) ,1)+diag(ones(1,N-1) ,-1);
Lap(1,2)=2;Lap(N,N-1)=2;

M1=delta*Lap/dx^2+eye(N)*(-1/dt -m);
M2=Lap/dx^2-eye(N);

maxn=[]; spread=[];

for idx=1:numel(tt)
    a=aa(idx);
    t=tt(idx);

    noise= randn(N,1)*sqrt(dt)*sigma0*sqrt(N)*0;
    wnext =(M2 -diag(n.^2))\(-a-noise/dt);
    nnext=M1\(-n/dt -n.^2.*w);
    n=nnext; w=wnext;

    maxn(end+1)=max(n);
    spread(end+1) = (max(n)-min(n))/mean(n);
    if t>tout
        tout=tout+10;
        subplot(2,1,1);
        plot(x,n,x,w);
        legend(' n' ,' w' ); xlabel( ' x' );
        title(sprintf( ' t=%g a=%g ' ,t, a));
        drawnow;
    end ;
end ;

subplot(2,1,2); hold on;
plot(aa,maxn);
at = interp1( spread , aa , 1); % revoir ca
d=delta;
ap = (3-2*sqrt(2-2*d))/(sqrt(3-2*sqrt(2-2*d)-d)*d);
disp("test")
disp(at)
disp(ap)
plot([at,at],[0,3], ' --b' );
xlabel(' a' ); ylabel( ' max(n)' );
title(sprintf( ' a_d=%g' ,at));