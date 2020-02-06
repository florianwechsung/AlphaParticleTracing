function [out] = kepler_inverse(x)
%KEPLER_INVERSE Summary of this function goes here
%   Detailed explanation goes here
if x>0.5
    s = 1-x;    
    s = 2*pi*s;
    s = (6*s)^(1/3);
else
    s = 2*pi*x;
    s = (6*s)^(1/3);
end
s2 = s^2;
out = s * (1 ...
        + s2*(1/60 ...
            + s2*(1/1400 ...
                + s2*(1/25200 ...
                    + s2*(43/17248000 ...
                        + s2*(1213/7207200000 ...
                        + s2*151439/12713500800000))))));
out = out/(2*pi);
if x>0.5
    out = 1-out;
end

