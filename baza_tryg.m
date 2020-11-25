function [ res ] = baza_tryg( x, K )
%Ortogonalna na przedziale [0, 2pi]
fik = @(x, k) 1/sqrt(pi)*cos((k/2).*x).*(mod(k,2) == 0) +       ...
              1/sqrt(pi)*sin(((k-1)/2).*x).*(mod(k-1,2) == 0) + ...
              (1/sqrt(2*pi)).*(k==1);   

res = zeros(K, 1);          
for i=1:1:K
   res(i) = fik(x, i);
end
          
end

