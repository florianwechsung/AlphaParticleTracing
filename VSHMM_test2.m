close all
eps = 1/1000;
rhs_slow_only = @(t, y) -eps * [y(1)*sin(y(2))*cos(y(2)); cos(y(2))^2];
rhs_full = @(t, y) rhs_slow_only(t, y) + [(y(1)*cos(y(2))+y(1)*sin(y(2))-(1/3)*y(1)^3*cos(y(2))^3)*cos(y(2)); -(cos(y(2)) + sin(y(2))-(1/3)*y(1)^2*cos(y(2))^3)*sin(y(2))]/eps;
y0 = [1/sin(0.25*pi);0.25*pi];


options = odeset('RelTol',1e-9, 'AbsTol', 1e-9);
[t, y] = ode45(rhs_full,[0 5000],y0, options);
slow = [y(:, 1).*sin(y(:, 2)) y(:, 1).*cos(y(:, 2))];
figure
hold on
% plot(t, y)
plot(t, slow)
%% 
tau = 0.05 * eps
delta = 0.01;%0.1 * tau/eps
delta_t = tau;
Delta_T = delta-delta_t;
iter = round(5000/(delta_t+Delta_T));
%iter = round(10/(delta_t));
[t, y] = FLAVORS(rhs_full,rhs_slow_only, y0, Delta_T, delta_t, iter);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 1, 0.1, 1e-4, 100);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 30, 1., 0.000294377, 10);
slow = [y(:, 1).*sin(y(:, 2)) y(:, 1).*cos(y(:, 2))];
%plot(t, y)
plot(t, slow)
