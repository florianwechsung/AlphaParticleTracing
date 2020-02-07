%p.19 fully nonlinear 5.2.2
close all
epsilon = 5*1e-4;
rhs_slow_only = @(t,y) [[[y(1);y(2)]/(sqrt(sqrt(y(1)^2 + y(2)^2))^3)+...
5*(y(1)/(sqrt(y(1)^2 + y(2)^2)) + y(3)/(sqrt(y(3)^2 + y(4)^2)))*[y(1);y(2)]/(sqrt(y(1)^2 + y(2)^2))];...
[[y(3);y(4)]/(sqrt(sqrt(y(3)^2 + y(4)^2)))^3 + y(3)*[y(3);y(4)]/(sqrt(y(3)^2 + y(4)^2))^2]];
  
rhs_full = @(t,y) rhs_slow_only(t, y) +[-y(2)*sqrt(y(1)^2 + y(2)^2)/epsilon; y(1)*sqrt(y(1)^2 + y(2)^2)/epsilon;...
    -y(4)*sqrt(y(3)^2 + y(4)^2)/(sqrt(2)*epsilon); y(3)*sqrt(y(3)^2 + y(4)^2)/(sqrt(2)*eps)];
y0 = [1;0;0;1];

options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
[t_dns, y_dns] = ode45(rhs_full,[0 3],y0, options);
slow_dns = sqrt(y_dns(:,1).^2 + y_dns(:,2).^2);

figure
hold on
plot(t_dns, slow_dns)

%5
tau = 0.05 * epsilon; 
%delta = 0.005
%iter = round(1/delta);
%[t, y] = FLAVORS(rhs_full,rhs_slow_only, y0, delta-tau, tau, iter);
Delta_T = 0.6; 
iter = round(3/Delta_T);
alpha = 10;
[t, y] = VSHMM(rhs_full,rhs_slow_only, y0, alpha, Delta_T, tau, iter);
% [t, y] = VSHMM(rhs_full,rhs_slow_only, y0, 30, 1., 0.000294377, 10);


slow = sqrt(y(:,1).^2 + y(:,2).^2);
plot(t, slow)