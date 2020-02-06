function [ts,ys] = FLAVORS(rhs_full,rhs_slow_only, y0, Delta_T, delta_t, niter)
%VSHMM Summary of this function goes here
%   Detailed explanation goes here

ts = [0];
ys = [y0'];

y = zeros(size(y0));
y(:) = y0(:);
t = 0;
counter = 0;
for i=1:niter
    tstep = delta_t;
    f1 = tstep*rhs_full(t, y);
    f2 = tstep*rhs_full(t+tstep/2, y+f1/2);
    f3 = tstep*rhs_full(t+tstep/2, y+f2/2);
    f4 = tstep*rhs_full(t+tstep, y+f3);
    y = y+(f1+2*(f2+f3)+f4)/6;
    
    t = t + tstep;
    sprintf('tstep=%f', tstep)
%     
    tstep = Delta_T;
    f1 = tstep*rhs_slow_only(t, y);
    f2 = tstep*rhs_slow_only(t+tstep/2, y+f1/2);
    f3 = tstep*rhs_slow_only(t+tstep/2, y+f2/2);
    f4 = tstep*rhs_slow_only(t+tstep, y+f3);
    y = y+(f1+2*(f2+f3)+f4)/6;
    t = t + tstep;
    sprintf('tstep=%f', tstep)


    sprintf('t=%f', t)
    counter = counter + 1;
    ts = [ts; t];
    ys = [ys; y'];
end

