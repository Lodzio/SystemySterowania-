% ZSS Projekt 1
clear all;
close all;

N = 10000; % Ilosc probek
gauss = 0 ; % Rodzaj wejscia
Her_ver = 6; % | 1 - podst | 2 - norm | 6 - tryg
u_przesuniecie = 1; % | 0 - U[0, 1] | 1 - U[0, 2pi]

% ---<<< Funkcja nieliniowa wspolczynniki
K_vec = [1, 2, 3]; 
K = length(K_vec); % Liczba Funkcji bazowych
% ---<<< Elementy odpowiedzi impulsowej
S = 5;                 
GammaV = linspace(5,1,S)';
 



% Wejscie
if gauss == 1
    U = randn(1,N); % Random input - N[0, 1]
else
    U = rand(1,N); % Random input - U[0, 1]
    if u_przesuniecie == 1
            U = (2*pi*rand(1,N)); % Random input - U[0, 2pi]
    end
end


% Sygnal Wk - niedostepny
Wk = fnlin(U, K_vec, Her_ver); 



Vk = zeros(1,N);

% Splot
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


% ---<<< Zaklocenie
Zk_variance = 0.1;
Zk_mean = 0;
Zk = Zk_variance.*randn(1,N) + Zk_mean;

%Yk = Vk + Zk;

Yk = Vk; % idealny przypadek bez zaklocenia



% === korelacja wzajemna ===
% Zdjecie skladowej stalej
U_centr = U - mean(U);
Y_centr = Yk - mean(Yk);

tau = S + 0; % ilos szukanych gamm
Gammy_est = zeros(tau, 1);
for i = 1:1:tau
    sum = 0;
    for j = 1:1:(N-i+1)
        sum = sum + U_centr(j)*Y_centr(j+i-1);
    end
    Gammy_est(i) = sum/(N-i+1);
end

if gauss
    disp('Wejscie Gaussowskie N');
else
    if u_przesuniecie == 1
        disp('Wejscie Jednostajne U[- \pi, \pi ]');
    else
        disp('Wejscie Jednostajne U[0, 1]');
    end
end
c = Gammy_est./GammaV;
str = ['Gammy estymowane / Gammy rzeczywiste Ratio'];
disp(str);
disp(c');

c = (GammaV / GammaV(1))' ;
str = ['Gammy / Gamma(1)'];
disp(str);
disp(c);

c = (Gammy_est / Gammy_est(1))' ;
str = ['Gammy_est/ Gamma_est(1)'];
disp(str);
disp(c);


L_nom = -500;
P_nom =  500;
Point = 1;
h_nom = 400;
A = zeros(1, K);
aloop = ones(1, K);
fexit = 1;
eps = 10^(-6);
for i=1:1:K
    fexit = 1;
    h = h_nom;
    L = L_nom;
    P = P_nom;
    loopiter = 1;
    while fexit
       Point = (L+P)/2;
       aloop(i) = Point + h;
       val_plus = Qfun(U, Yk, Gammy_est, aloop, Her_ver);
       aloop(i) = Point - h;
       val_minus = Qfun(U, Yk, Gammy_est, aloop, Her_ver);
       if val_plus >= val_minus
           P = Point + h;
       else
           L = Point - h; 
       end
       if (P - L) <= eps
           fexit = 0;
           str = ['Obliczono A',num2str(i), ' = ',num2str(aloop(i)),' po ', num2str(loopiter) , ...
                  ' iteracjach, Q = ', num2str(Qfun(U, Yk, Gammy_est, aloop, Her_ver))];
           disp(str);
       end
       loopiter = loopiter +1;
       h = h*0.5;
    end
    A(i) =P;
    
end
