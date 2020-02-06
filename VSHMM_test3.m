close all
eps = 1/3400;
rhs_slow_only = @(t, y) 0.25 * y + 5 * y(1) * y / sqrt(y(1)^2 + y(2)^2);
rhs_full = @(t, y) rhs_slow_only(t, y) + [-y(2); y(1)]/eps;
y0 = [1;0];

options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
[t, y] = ode45(rhs_full,[0 5],y0, options);
slow_dns = sqrt(y(:,1).^2 + y(:,2).^2);
%%
% figure
hold on
% plot(t, y)
plot(t, slow_dns)
%% 
tau = 0.05 * eps

iter = round(1/delta);
%[t, y] = FLAVORS(rhs_full,rhs_slow_only, y0, delta-tau, tau, iter);
Delta_T = 1.0;
iter = round(5/Delta_T);
alpha = 50;
[t, y] = VSHMM(rhs_full,rhs_slow_only, y0, alpha, Delta_T, tau, iter);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 30, 1., 0.000294377, 10);

slow = sqrt(y(:,1).^2 + y(:,2).^2);
%%
plot(t, slow)