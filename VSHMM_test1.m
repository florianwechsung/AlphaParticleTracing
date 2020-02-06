close all
eps = 1/1000;
% rhs_slow_only = @(t, y) 0.25 * y + 5 * y(1) * y / sqrt(y(1)^2 + y(2)^2);
% rhs_full = @(t, y) rhs_slow_only(t, y) + [-y(2); y(1)]/eps;
% y0 = [1;0];

rhs_slow_only = @(t, y) [-eps*y(2); 0];
rhs_full = @(t, y) rhs_slow_only(t, y) + [0; y(1) + y(2) - (1./3) * y(2)^3]/eps;
y0 = [1;1];


options = odeset('RelTol',1e-9, 'AbsTol', 1e-9);
[t, y] = ode45(rhs_full,[0 5000],y0, options);
slow = sqrt(y(:,1).^2 + y(:,2).^2);
figure
hold on
plot(t, y)
% 
tau = 0.05 * eps
delta = 0.01;
delta_t = tau;
Delta_T = delta-delta_t;
iter = round(5000/(delta_t+Delta_T));
[t, y] = FLAVORS(rhs_full,rhs_slow_only, y0, Delta_T, delta_t, iter);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 1, 0.1, 1e-4, 100);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 30, 1., 0.000294377, 10);
plot(t, y)