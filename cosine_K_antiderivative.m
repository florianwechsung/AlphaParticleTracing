function [K_integral] = cosine_K_antiderivative(t,Delta_T)
%COSINE_K_ANTIDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here
K_integral = Delta_T * (t/Delta_T - (1/(2*pi)) * sin(2*pi*(t/Delta_T)));
end

