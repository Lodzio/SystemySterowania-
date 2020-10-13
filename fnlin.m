function y = fnlin(u,k_vec)
% function return vector of the output which is dot product of
% hermite polynomians and vector k_vec.
% u is the input vector of the function
% k_vec is the vector of multipliers for base functions
% The length of vector y is equal to length of vector u

k = length(k_vec);
n = length(u);
y = zeros(1,n);

for i = 1:n
   hermite_v = hermite(u(i), k-1);
   y(i) = k_vec * hermite_v;
end
end
