clear all;
close all;

u = [2, 2];

min1 = getMin(@(x) Q([x 10]), -10, 10);
min2 = getMin(@(x) Q([min1 x]), -10, 10);
Pmin = [min1, min2]
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


function y = Q(u)
    y=(u(1)-4)^2 + (u(2)-4)^2;
end