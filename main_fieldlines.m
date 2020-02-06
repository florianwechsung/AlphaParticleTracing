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
%                       Field line tracing parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=5;
N0 = 250;%N0 is the number of toroidal turns - 500 default number so far
Nstep = 40;%Nstep is the number of Runge-Kutta steps per turn
N1 = Nstep*N0;%N1 is the number of steps for Runge-Kutta
Ngrid=80;%Ngrid is the number of grid points for the x and y regular grid on which one interpolates
t0 = 0;%Starting toroidal angle 
tf = N0*2*pi;%Ending toroidal angle
y01 = [1.136 1.183 1.242 1.285 1+epsilon];%Array of starting x
y02 = 0;%Starting y

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
%                   Field line tracing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize all the quantities
 U=zeros(N1+1,b);
 V=zeros(N1+1,b);

for i=1:b
[t,y] = ode_RK4(@MyFieldLines, [t0 tf], [y01(i) y02], N1,C,A,D(1),D(2),D(3),Btin,Rin);
 U(:,i)=y(:,1);
 V(:,i)=y(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   Plot results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
[X,Y] = meshgrid(0:.001:1+epsilon+0.5,-kappa*epsilon-0.1:.001:kappa*epsilon+0.1);

Z = C*X.^4/8+A*(X.^2.*log(X)/2-X.^4/8)+...
    D(1)+D(2)*X.^2+D(3)*(X.^4-4*X.^2.*Y.^2);
I = find(Z>0);
Z(I)=NaN;

figure(1)
pcolor(X, Y, Z);
shading interp
contourf(X,Y,Z);
hold on
plot(U(1:Nstep:N1-Nstep+1,1),V(1:Nstep:N1-Nstep+1,1),'Color','black','Marker','.','LineStyle', 'none')
hold on
plot(U(1:Nstep:N1-Nstep+1,2),V(1:Nstep:N1-Nstep+1,2),'Color','black','Marker','.','LineStyle', 'none')
hold on
plot(U(1:Nstep:N1-Nstep+1,3),V(1:Nstep:N1-Nstep+1,3),'Color','black','Marker','.','LineStyle', 'none')
hold on
plot(U(1:Nstep:N1-Nstep+1,4),V(1:Nstep:N1-Nstep+1,4),'Color','black','Marker','.','LineStyle', 'none')
hold on
plot(U(1:Nstep:N1-Nstep+1,5),V(1:Nstep:N1-Nstep+1,5),'Color','black','Marker','.','LineStyle', 'none')
axis equal
grid on

