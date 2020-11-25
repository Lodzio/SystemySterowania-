function y = Qfun(U, Y, est_G, est_A, hv)
% function return vector of the output which is dot product of
% hermite polynomians and vector k_vec.
% u is the input vector of the function
% k_vec is the vector of multipliers for base functions
% The length of vector y is equal to length of vector u



N = length(U); % Number of samples
Wk = fnlin(U, est_A, hv); % Calculating Wk Vector

S = length(est_G);                 % Amount of Gamma's
GammaV = est_G;

Vk = zeros(1,N);

% Convolution calculation
for i=1:N
    phi = zeros(1, S);
    for j=1:S
        if i - j + 1 <= 0    % if idx less t
            skl_mi = 0;
        else
            skl_mi = Wk(i-j+1);
        end  
        phi(j) = skl_mi;
    end
   Vk(i) = phi*GammaV;
end
y = ((Y-Vk)*(Y-Vk)')/N;
