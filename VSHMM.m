function [ts,ys] = VSHMM(rhs_full,rhs_slow_only, y0, alpha, Delta_T, delta_t, niter)
%VSHMM Summary of this function goes here
%   Detailed explanation goes here

ts = [0];
ys = [y0'];

y = zeros(size(y0));
y(:) = y0(:);
t = 0;
counter = 0;
for i=1:niter
    tlocal = 0.;
    while tlocal < Delta_T
        tstep = min(delta_t, Delta_T-tlocal);
        %tstep = delta_t;
        %y = y + tstep * rhs_full(t, y);
        f1 = tstep*rhs_full(t, y);
        f2 = tstep*rhs_full(t+tstep/2, y+f1/2);
        f3 = tstep*rhs_full(t+tstep/2, y+f2/2);
        f4 = tstep*rhs_full(t+tstep, y+f3);
        y = y+(f1+2*(f2+f3)+f4)/6;
    
        tlocal = tlocal + tstep;
        t = t + tstep;
        sprintf('tstep=%f', tstep)
        
        
        h_t = alpha * delta_t * cosine_K(cosine_K_antiderivative_inverse(mod(t, Delta_T), Delta_T), Delta_T);
        tstep = min(h_t, Delta_T-tlocal);
        %tstep = h_t;
%         y = y + tstep * rhs_slow_only(t, y);
        f1 = tstep*rhs_slow_only(t, y);
        f2 = tstep*rhs_slow_only(t+tstep/2, y+f1/2);
        f3 = tstep*rhs_slow_only(t+tstep/2, y+f2/2);
        f4 = tstep*rhs_slow_only(t+tstep, y+f3);
        y = y+(f1+2*(f2+f3)+f4)/6;
        tlocal = tlocal + tstep;
        t = t + tstep;
        sprintf('tstep=%f', tstep)
        sprintf('t=%f', t)
        counter = counter + 1;
        ts = [ts; t];
        ys = [ys; y'];
    end
    counter
end

