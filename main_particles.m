clear

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
T_particleTracing = 2000*T;%Total simulation time: 2000 cyclotron periods (just for checking trajectories at the moment)
dT = pi/(32*omega_c);%Size of the time step for numerical ode solver
MM = T_particleTracing/dT;%Number of time steps
 
[tpart,ypart] = ode_RK4(@MyParticleTrajectory, [0 T_particleTracing], [1+epsilon/2 5e5 0 1e5 0 0], MM,C,A,D(1),D(2),D(3),Btin,Rin,omega_c);%Initial conditions for nice banana orbit
%[tpart,ypart] = ode_RK4(@MyParticleTrajectory, [0 T_particleTracing], [1+epsilon/2 1e3 0 1e5 0 0], MM,C,A,D(1),D(2),D(3),Btin,Rin,omega_c);%Initial conditions for passing trajectory
Rpart=ypart(:,1);%Radial coordinate in cylindrical coordinate system
phipart=ypart(:,3);%Angle coordinate in cylindrical coordinate system
Zpart=ypart(:,5);%Z coordinate in cylindrical coordinate system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Plot results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct equilibrium flux surfaces
[X,Y] = meshgrid(0:.001:1+epsilon+0.5,-kappa*epsilon-0.1:.001:kappa*epsilon+0.1);

Z = C*X.^4/8+A*(X.^2.*log(X)/2-X.^4/8)+...
    D(1)+D(2)*X.^2+D(3)*(X.^4-4*X.^2.*Y.^2);
I = find(Z>0);
Z(I)=NaN;

% Plot flux surfaces and trajectories

figure(1)
pcolor(X, Y, Z);
shading interp
contourf(X,Y,Z);
hold on
plot(Rpart,Zpart,'Color','red','Marker','o','MarkerSize',2,'LineStyle', 'none')
axis equal
grid on

