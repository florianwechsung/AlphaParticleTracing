function [K] = cosine_K(t, Delta_T)
%COSINE_K Summary of this function goes here
%   Detailed explanation goes here
K = (1/Delta_T) * (1 + cos(2*pi*(t/Delta_T-0.5)));
end

