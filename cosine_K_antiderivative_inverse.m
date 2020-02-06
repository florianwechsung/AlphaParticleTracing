function [t] = cosine_K_antiderivative_inverse(K,Delta_T, options)
%COSINE_K_ANTIDERIVATIVE_INVERSE Summary of this function goes here
%   Detailed explanation goes here
residual = @(t) deal(K-cosine_K_antiderivative(t, Delta_T), -cosine_K(t, Delta_T));
t = kepler_inverse(K/Delta_T)*Delta_T;
if abs(K-cosine_K_antiderivative(t, Delta_T)) > 1e-4*Delta_T
    t = fsolve(residual, t, options);
end
end