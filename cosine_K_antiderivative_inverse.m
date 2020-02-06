function [t] = cosine_K_antiderivative_inverse(K,Delta_T)
%COSINE_K_ANTIDERIVATIVE_INVERSE Summary of this function goes here
%   Detailed explanation goes here
residual = @(t) deal(K-cosine_K_antiderivative(t, Delta_T), -cosine_K(t, Delta_T));
options = optimset('Jacobian', 'on', 'Display','off');
t = fsolve(residual, K, options);
end