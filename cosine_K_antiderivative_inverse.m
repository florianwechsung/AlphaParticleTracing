function [t] = cosine_K_antiderivative_inverse(K,Delta_T)
%COSINE_K_ANTIDERIVATIVE_INVERSE Summary of this function goes here
%   Detailed explanation goes here
residual = @(t) K-cosine_K_antiderivative(t, Delta_T);
options = optimset('Display','off');
t = fsolve(residual, K, options);
end