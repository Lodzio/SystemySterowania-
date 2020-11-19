% ZSS Projekt 1
clear all;
close all;

N = 5000; % Number of samples
U = rand(1,N); % Random input - u[0, 1]
%K_vec = [1, 2, 3, 4, 5]; % Vector of multipliers of base function
                         % Must be size [1 x K] because of dot product
                         % in nonlinear function

 K_vec = [5, 2]; % mniejszy przypadek zeby zobaczyc czy dziala

K = length(K_vec); % Number of Base functions

Wk = fnlin(U, K_vec); % Calculating Wk Vector

% ---<<< Impulse response (Gamma's)
S = 5;                 % Amount of Gamma's
numerator = 1;          %   1
denumerator = [1, 10];  % s + 10
h = tf(numerator, denumerator);
t = linspace(0,1,S);
[im_y, im_x] = impulse(h, t);
%GammaV = im_y;
GammaV = linspace(5,1,S)';

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



% korelacja wzajemna
U_centr = U - mean(U);
Y_centr = Yk - mean(Yk);

tau = S + 0; % ilos szukanych gamm
Gammy_est = zeros(tau,1);
%{
for i = 1:1:tau
    sum = 0;
    for j = 1:1:(N-i+1)
        sum = sum + (U_centr(j)*Y_centr(j+i-1));
    end
    Gammy_est(i) = sum/(N-i+1);
end
%}
for i = 1:tau
    sum = U_centr(1:N-i+1)*Y_centr(i:end)';
    Gammy_est(i) = sum/(N-i+1);
end
Gammy_est ./ GammaV
Gammy_est = GammaV;
%Gammy_est = GammaV*-2.6;

% Gammy_est / Gammy_est(1)
% (GammaV / GammaV(1))'

% a = Qfun(U, Yk, Gammy_est', A)

L_nom = -10;
P_nom =  10;
goldenRatio = (( sqrt(5)-1 ) / 2);
Point = 1;
A = zeros(1, K);
aloop = ones(1,K);
fexit = 1;
eps = 10^(-2);
for ind=1:100
    i = rem(ind, K)+1
    fexit = 1;
    L = L_nom;
    P = P_nom;
    loopiter = 0;
    while fexit
       h = (P-L)*goldenRatio;
       aloop(i) = L + h;
       val_plus = Qfun(U, Yk, Gammy_est, aloop);
       aloop(i) = P - h;
       val_minus = Qfun(U, Yk, Gammy_est, aloop);
       if val_plus >= val_minus
           P = L + h;
       else
           L = P - h; 
       end
       if (abs(P - L)) <= eps
           fexit = 0;
           str = ['Obliczono A',num2str(i), ' = ',num2str(aloop(i)),' po ', num2str(loopiter) , ...
                  ' iteracjach, Q = ', num2str(Qfun(U, Yk, Gammy_est, aloop))];
           disp(str);
       end
       loopiter = loopiter +1;
    end
    A(i) = aloop(i);
    
end

A./K_vec



%{
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
Teta = (Phi'*Phi)^-1*Phi'*Yk';
M = reshape(Teta, [K, S]);
[P, D, Q] = svd(M);
K_hat = P(:,1);
sigma1 = D(1);
G_hat = Q(:,1); % bo svd matlaba zwraca PDQ'

% ratio
(GammaV(2:end)./GammaV(1:end-1))'
(G_hat(2:end)./G_hat(1:end-1))'

K_vec(2:end)./K_vec(1:end-1)
(K_hat(2:end)./K_hat(1:end-1))'
%}