clear all;
close all;


% policzenie optymalnych u1, u2
sterowanie = getMiniumum_u1u2(1)

% suma kwadratow sterowan
sterowanie*sterowanie'

% policzenie wyjœcia systemu dla takich sterowañ
wyjscie = system(sterowanie')

% funkcja celu wartosc dla policzonego sterowania
Q_basic(wyjscie)


function q = Q_basic(y)
    % u = [ u1, u2]
    q = (y(1) - 4).^2 + (y(2) - 4).^2;
end

%{

Rysowanie 3d do zrobienia

function drawQfun(r, N)
    % promien zmiennosci r1,r2 tzn r1 w [-r1, r1] etc.
    figure();
    y1vec = linspace(-r(1), r(1), N);
    y2vec = linspace(-r(2), r(2), N);
    Y = [y1vec', y2vec'];
    res =  Q_basic(Y);
    surf ( y1vec, y2vec, res);
    xlabel("x");
    ylabel("t");
    zlabel("T");
    %title({"k = " + num2str(k) + " | T1 = " + num2str(T1) + " | T2 = " + num2str(T2)},'FontSize',20);  

end
%}

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
 
 
 
function min = getMin(fun, Lstart, Pstart)
    stopValue = 1e-7;
    E = 1;
    a = Lstart;
    b = Pstart;
    while E > stopValue
        center = (a+b)/2;
        min = center;
        P = center + E;
        L = center - E;
        Pvalue = fun(P);
        Lvalue = fun(L);
        if Pvalue > Lvalue
            b = P;
        else
            a = L;
        end
        E = E/2;
    end
end

function q = Q(u)
    A=[0.5 0; 0 0.25];
    B=[1 0; 0 1];
    H=[0 1; 1 0];
    K = (eye(size(A))-A*H)^-1*B;
    y = K*u.';
    q=(y(1)-4)^2 + (y(2)-4)^2;
end
