close all
clear

%---------------- Parameter und lineares GLG von Antoine übernommen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Parameters defining equilibrium
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon=0.32; % Inverse aspect ratio of equilibrium
kappa=1.7; % Elongation of equilibrium
delta=0.33; % Triangularity of equilibrium
A = -0.2; % A parameter for Solov'ev profile for F
C=1-A; % C parameter for Solov'ev profile for p
Rout = 1+epsilon; % R at the outboard midplane
Rin = 1-epsilon; %R at the inboard midplane
Rtop = 1-delta*epsilon; % R at top point
Ztop = kappa*epsilon; % Z at top point
Btin = 1.5; % Toroidal magnetic field at inboard midplane; choose 1 for large change in the poloidal angle after one turn; choose 60 for small change in the poloidal angle after one turn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Construct and solve linear system for free constants c1,
%               c2, c3, called D(1), D(2), and D(3) here
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Amat is the matrix corresponding to imposing the constraint conditions on
% the homogeneous solutions

Amat = [1 Rout^2 Rout^4
    1 Rin^2 Rin^4
    1 Rtop^2 Rtop^4-4*Rtop^2*Ztop^2];

% B is the matrix corresponding to imposing the constraint conditions on
% the inhomogeneous terms

B= -[C/8*Rout^4+A*(Rout^2*log(Rout)/2-Rout^4/8)
     C/8*Rin^4+A*(Rin^2*log(Rin)/2-Rin^4/8)
     C/8*Rtop^4+A*(Rtop^2*log(Rtop)/2-Rtop^4/8)];

% Solve for free coefficients 
 
D=Amat\B;

c1 = D(1);
c2 = D(2);
c3 = D(3);
Bt = Btin;
xb = Rin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Particle trajectories
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = 2;%Particle charge, in units of e
m = 6.64e-27;%Particle mass
omega_c = q*1.6e-19*Btin/m;%Cyclotron angular frequency at the inboard midplane
T = 2*pi/omega_c;%Cyclotron period
M = 1e7;%Approximate number of cyclotron periods followed for the particle trajectories
T_particleTracing = 2*T;%Total simulation time: 2000 cyclotron periods (just for checking trajectories at the moment)
dT = pi/(32*omega_c);%Size of the time step for numerical ode solver
MM = T_particleTracing/dT;%Number of time steps

%------------------------ Ende Antoine übernommen

rhs_slow_only = @(t, y) [y(2)
    y(1)*y(4)^2
    y(4)
    -2*y(2)/y(1)*y(4)
    y(6)
    0];

rhs_full = @(t, y) [y(2)
    omega_c*(y(1)*y(4)*Bz(y(1),y(5),A,C,c1,c2,c3)...
    -y(6)*Bphi(y(1),y(5),A,C,c1,c2,c3,Bt,xb))+y(1)*y(4)^2
    y(4)
    omega_c/y(1)*(y(6)*Br(y(1),y(5),A,C,c1,c2,c3)...
    -y(2)*Bz(y(1),y(5),A,C,c1,c2,c3))-2*y(2)/y(1)*y(4)
    y(6)
    omega_c*(y(2)*Bphi(y(1),y(5),A,C,c1,c2,c3,Bt,xb)...
    -y(1)*y(4)*Br(y(1),y(5),A,C,c1,c2,c3))];
y0 = [1+epsilon/2 5e5 0 1e5 0 0]'; % for banana orbit
%y0 = [1+epsilon/2 1e3 0 1e5 0 0]'; % for passing trajectory

options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
[t_dns, y_dns] = ode45(rhs_full,[0 T_particleTracing],y0, options);
slow_dns = sqrt(y_dns(:,1).^2 + y_dns(:,2).^2);

%%
 figure
hold on
% plot(t, y)
plot(t_dns, slow_dns)

%% 
%tau = 0.05 * eps
tau = dT;
% delta = 0.01
% iter = round(1/delta);
%[t, y] = FLAVORS(rhs_full,rhs_slow_only, y0, delta-tau, tau, iter);
Delta_T = T_particleTracing/5;
iter = round(T_particleTracing/Delta_T);
alpha = 10;
[t, y] = VSHMM(rhs_full,rhs_slow_only, y0, alpha, Delta_T, tau, iter);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 30, 1., 0.000294377, 10);

slow = sqrt(y(:,1).^2 + y(:,2).^2);
%%
plot(t, slow)