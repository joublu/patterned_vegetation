% Turing instability of the modified Klausmeier model (Lecture 24, no advection term)

clear all
close all
clc

% modified Klausmeier (only diffusion)

% parameters trees
% a: 0.077 - 0.23
% m=0.045;

% parameters grass
% a: 0.94 - 2.81
m=0.45; % plant mortality (natural and grazing)
% aC=0.9; % saddle-node bifurcation point
d=0.04; % diffusion coefficient for plants


figure()
hold on
box on
plot([0 10],[0 0],'k:')
kk=0:0.1:10; % possible values for the eigenvalues k^2

% Determinant of the characteristic matrix for decreasing a
for a=[1 0.99 0.9857 0.98 0.97 0.96 0.95 0.94] %grass
% Homogeneous state (coexistence)

% E minus (stable when it exists, a>2m)
ns=a./(2*m)+sqrt((a/(2*m)).^2-1);
ws=m./ns;

% Jacobian
J11=m;
J12=ns.^2;
J21=-2*m;
J22=-1-ns.^2;

TrJs=J11+J22;
DetJs=J11*J22-J21*J12;

disp(TrJs)
disp(DetJs)

% Determinant of the charachteristic matrix Mk
DetMk=@(kk) d*kk.^2-(J11+d*J22)*kk+DetJs;

plot(kk,DetMk(kk))
axis([0 10 -0.4 1])
xlabel('k^2')
ylabel('Det(M_k)')

caption = sprintf('Det(M_k) for a= %.3f', a);
title(caption, 'FontSize', 14);
pause

end

%% Pattern formation
% now we fix a domain and we should check if at least one eigenvalue falls
% into the TIR

close all
% 1D domain, hNbc
L=10; %domain length
a=0.94;

% E minus (stable when it exists, a>2m)
ns=a./(2*m)+sqrt((a/(2*m)).^2-1);
ws=m./ns;

% Jacobian
J11=m;
J12=ns.^2;
J21=-2*m;
J22=-1-ns.^2;

DetJs=J11*J22-J21*J12;

% Determinant of the charachteristic matrix Mk
DetMk=@(kk) d*kk.^2-(J11+d*J22)*kk+DetJs;

figure()
hold on
box on
plot([0 10],[0 0],'k:')
plot(kk,DetMk(kk))
axis([0 10 -0.4 1])
xlabel('k^2')
ylabel('Det(M_k)')

caption = sprintf('Det(M_k) for a= %.3f', a);
title(caption, 'FontSize', 14);

% we look at the position of some eigenvalues
nn=1:30;
kk2_L=(pi*nn/L).^2;
plot(kk2_L,0*kk2_L,'*')
pause

% eigenfunction associated to the most unstable mode
% we select nn from the figure
[~,mum_i]=min(DetMk(kk2_L));
nn_mu=nn(mum_i);
xx=0:0.01:L;
figure()
plot(xx,cos(pi*nn_mu*xx/L)) % the pattern forming will look like the eigenfunction associated to the most unstable mode (at least close to the bifurcation!)

%%
% estimate of aT (critical value for Turing instability)
ns_a=@(a) a./(2*m)+sqrt((a/(2*m)).^2-1);
J11_a=m;
J12_a=@(a) ns_a(a).^2;
J21_a=-2*m;
J22_a=@(a) -1-ns_a(a).^2;

%TrJs_a=@(a) J11_a+J22_a(a);
DetJs=@(a) J11_a*J22_a(a)-J21_a*J12_a(a);
discriminant=@(a) (J11_a+d*J22_a(a)).^2-4*d*DetJs(a);

% aT(grass)=0.9857
aT=fzero(discriminant,1);
disp(aT)

%% Bifurcation diagram 
close all
% Desert state (always stable)
E0_n=0;
%E0_w=a; 

% E plus (unstable when it exists)
Ep_n=@(a) a./(2*m)-sqrt((a/(2*m)).^2-1);
%Ep_w=m./Ep_n;

% E minus (stable when it exists, a<2m)
Em_n=@(a) a./(2*m)+sqrt((a/(2*m)).^2-1);
%Em_w=m./Em_n;

figure()
hold on
box on
a=2*m:0.001:1.5;
plot([0 a(end)],[0,0],'k','LineWidth',2)
plot(a,Ep_n(a),'b:','LineWidth',2)
plot(a,Em_n(a),'b','LineWidth',2)
xlabel('a')
ylabel('n_{eq}')
axis([0.8,1.2,0,2.5])
plot([2*m 2*m],[0 1],'k:')%aC
plot([aT aT],[0 5],'k:')%aC

pause
a_pattern=2*m:0.001:aT;
plot(a_pattern,Em_n(a_pattern),'r:','LineWidth',2)
a_nopattern=aT:0.001:1.5;
plot(a_nopattern,Em_n(a_nopattern),'r','LineWidth',2)

legend('E_{0}','E_{+}','E_{-}','','','')