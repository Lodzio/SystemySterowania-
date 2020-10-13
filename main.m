% ZSS Projekt 1
clear all;
close all;

N = 1000; % Number of samples
U = rand(1,N); % Random input - u[0, 1]
K_vec = [1, 2, 3, 4, 5]; % Vector of multipliers of base function
                         % Must be size [1 x K] because of dot product
                         % in nonlinear function

K = length(K_vec); % Number of Base functions

Wk = fnlin(U, K_vec); % Calculating Wk Vector

% ---<<< Impulse response (Gamma's)
S = 11;                 % Amount of Gamma's
numerator = 1;          %   1
denumerator = [1, 10];  % s + 10
h = tf(numerator, denumerator);
t = linspace(0,1,S);
[im_y, im_x] = impulse(h, t);
GammaV = im_y;

Vk = zeros(1,N);

% Convolution calculation
for i=1:N
    sumi = 0;
    for j=1:S
        if i - j + 1 <= 0    % if idx less t
            skl_mi = 0;
        else
            skl_mi = Wk(i-j+1);
        end  
        sumi = sumi + (GammaV(j)*skl_mi);
    end
   Vk(i) = sumi;
end

% ---<<< Noise 
Zk_variance = 0.1;
Zk_mean = 0;
Zk = Zk_variance.*randn(1,N) + Zk_mean;

Yk = Vk + Zk;