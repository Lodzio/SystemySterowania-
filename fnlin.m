function y = fnlin(u,k_vec, hv)
% function return vector of the output which is dot product of
% hermite polynomians and vector k_vec.
% u is the input vector of the function
% k_vec is the vector of multipliers for base functions
% The length of vector y is equal to length of vector u

k = length(k_vec);
n = length(u);
y = zeros(1,n);


if hv <= 5
    for i = 1:n
       hermite_v = hermite(u(i), k-1, hv);
       y(i) = k_vec * hermite_v;
    end
else
    for i = 1:n
        tmpv = baza_tryg(u(i), k);
        y(i) = k_vec*tmpv;
    end  
end
    
