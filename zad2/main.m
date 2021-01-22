clear all;
close all;
r = 1;

u1min = getMin(@(u1) getQMinForU1(u1, r), -r, r);
getMinU2Comparator = @(u2) comparator([u1min; u2]);
u2min = getMin(getMinU2Comparator, -sqrt(r^2 - u1min^2), sqrt(r^2 - u1min^2));

uMin = [u1min, u2min]
% policzenie optymalnych u1, u2
sterowanie = getMiniumum_u1u2(1)

% rysowanie
drawQfun([10 5], 51);

% suma kwadratow sterowan
% sterowanie*sterowanie'

% policzenie wyj?cia systemu dla takich sterowa?
% wyjscie = system(sterowanie')

% funkcja celu wartosc dla policzonego sterowania
% Q_basic(wyjscie)


function q = Q_basic(y)
    q = (y(1) - 4).^2 + (y(2) - 4).^2;
end

function drawQfun(r, N)
    % promien zmiennosci r1,r2   
    y1vec = linspace(-r(1), r(1), N);
    y2vec = linspace(-r(1), r(1), N);
    u1vec = linspace(-r(2), r(2), N);
    u2vec = linspace(-r(2), r(2), N);
    
    ValueMY = zeros(N, N);
    ValueMU = zeros(N, N);

    for i=1:1:N
       for j=1:1:N
           ValueMY(i, j) = Q_basic([y1vec(i) y2vec(j)]);
           ValueMU(i, j) = comparator([u1vec(i); u2vec(j)]);
       end
    end
    figure();
    surf ( y1vec, y2vec, ValueMY);
    xlabel("y1");
    ylabel("y2");
    zlabel("Q(y1, y2)");
    figure();
    hold on;
    surf ( u1vec, u2vec, ValueMU);
    xlabel("u1");
    ylabel("u2");
    zlabel("Q(u1, u2)");
    
    R = 1;
    H = 54;
    [X,Y,Z] = cylinder(R);
    Z = Z*H;
    mesh(X,Y,Z,'FaceAlpha','0.5');
    hold off;
    grid on;

end

function y = system(u)
    % na podstawie [u1, u2] liczymy wyjscia sys [y1, y2]
    A=[0.5 0; 0 0.25];
    B=[1 0; 0 1];
    H=[0 1; 1 0];
    K = pinv(eye(2)-A*H)*B;
    y = K*u;
end

function res = getMinimum_u2(u1, ru2)
    % u1 - ustalone poziom nadrzedny
    % ru2 - promien po ktorym idziemy ze zmiennoscia u2 in [-ru2, ru2]
    % res - optymalne u2 przy zadanym u1 z dokladnoscia eps
    
    eps = 1e-7;
    L = -ru2;
    R = ru2;
    center = (L + R)/2;
    E1 = ru2/2;
    while E1 > eps
        center = (L + R)/2;
        Lcenter = center - E1;
        Rcenter = center + E1;
        
        LcenterY = system([u1; Lcenter]); % policzenie wyjsc systemu
        RcenterY = system([u1; Rcenter]); % policzenie wyjsc systemu
        
        QL = Q_basic(LcenterY'); % kryterium jakosci
        QR = Q_basic(RcenterY');
        if(QL >= QR)
            L = Lcenter;
        else
            R = Rcenter;
        end
        E1 = E1/2;
    end
    res = center;
end

function res = getMiniumum_u1u2(r)
    % r = promien kola
    % ograniczenie postaci u1^2 + u2^2 < r^2
    Lu1 = -r;
    Ru1 =  r;
    eps = 1e-7;
    E1 = r/2;
    while E1 > eps
        center_u1 = (Lu1+Ru1)/2;
        Lcenter_u1 = center_u1 - E1;
        Rcenter_u1 = center_u1 + E1;
        
        % optymalizacja po u2
        L_ru2 = sqrt(r^2 - Lcenter_u1^2); % promien zmiennosci u2
        R_ru2 = sqrt(r^2 - Rcenter_u1^2);
        
        % Policzenie najlepszego mozliwego u2 przy ustalonym u1
        L_u2_min = getMinimum_u2(Lcenter_u1, L_ru2);
        R_u2_min = getMinimum_u2(Rcenter_u1, R_ru2);
        
        %Policzenie wyjscia systemu dla zadanych u1 u2
        L_Y = system([Lcenter_u1; L_u2_min]);
        R_Y = system([Rcenter_u1; R_u2_min]);
        
        QL = Q_basic(L_Y');
        QR = Q_basic(R_Y');
        if(QL >= QR)
            Lu1 = Lcenter_u1;
        else
            Ru1 = Rcenter_u1;
        end
        E1 = E1/2; 
        res = [center_u1 R_u2_min];
    end
end

 % --- END Karol 
 
 
 
function min = getMin(comparator, Lstart, Pstart)
    stopValue = 1e-7;
    a = Lstart;
    b = Pstart;
    E = Pstart/2;
    while E > stopValue
        center = (a+b)/2;
        P = center + E;
        L = center - E;
        Pvalue = comparator(P);
        Lvalue = comparator(L);
        if Pvalue >= Lvalue
            b = P;
        else
            a = L;
        end
        E = E/2;
    end
    min = (a+b)/2;
end

function q = comparator(u)
    y = system(u);
    q = Q_basic(y);
end

function q = getQMinForU1(u1, r)
    Lstart = -sqrt(r^2 - u1^2);
    Pstart =  sqrt(r^2 - u1^2);
    u2 = getMin(@(u2) comparator([u1; u2]), Lstart, Pstart);
    q = comparator([u1; u2]);
end