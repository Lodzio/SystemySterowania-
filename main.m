% ZSS Projekt 1
clear all;
close all;

N = 10000; % Number of samples
U = rand(1,N); % Random input - u[0, 1]
K_vec = [1, 2, 3, 4, 5]; % Vector of multipliers of base function
                         % Must be size [1 x K] because of dot product
                         % in nonlinear function

K_vec = [1, 2, 3]; % mniejszy przypadek zeby zobaczyc czy dziala

K = length(K_vec); % Number of Base functions

Wk = fnlin(U, K_vec); % Calculating Wk Vector

% ---<<< Impulse response (Gamma's)
S = 6;                 % Amount of Gamma's
numerator = 1;          %   1
denumerator = [1, 10];  % s + 10
h = tf(numerator, denumerator);
t = linspace(0,1,S);
[im_y, im_x] = impulse(h, t);
GammaV = im_y;

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

% ---<<< Noise 
Zk_variance = 0.1;
Zk_mean = 0;
Zk = Zk_variance.*randn(1,N) + Zk_mean;

%Yk = Vk + Zk;
Yk = Vk; % idealny przypadek bez zaklocenia

% parametric identification
Phi = [];
for i = 1:N
    phi = [];
    for j = 0:S-1
        if i - j > 0
            h = hermite(U(i-j), K - 1);
        else
            h = zeros(K, 1);
        end
        phi=[phi, h.'];
    end
   Phi = [Phi; phi]; % bo rozmiar ma byc N x t*(s-1)
end
Teta = (Phi'*Phi)^-1*Phi'*Yk.';
M = reshape(Teta, [K, S]);
[P, D, Q] = svd(M);
K_hat = P(:,1);
sigma1 = D(1);
G_hat = Q(:,1); % bo svd matlaba zwraca PDQ'

% ratio
(GammaV(2:end)./GammaV(1:end-1))';
(G_hat(2:end)./G_hat(1:end-1))';

K_vec(2:end)./K_vec(1:end-1);
(K_hat(2:end)./K_hat(1:end-1))';

%%%%% zad 3
phi = zeros(N, length(K_vec))
for i = 1:length(U)
    for j = 1:length(K_vec)
        if i-j+1 > 0
            phi(i, j) = U(i-j+1);
        else 
            phi(i, j) = 0;
        end
    end
end

gamma = ((phi'*phi)^-1)*phi'*Yk'
K_vec./gamma'