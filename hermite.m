function y = hermite(x,k, her_version)
% Return vector with K order hermite function values at point x
% size of the vector is [k x 1]

if her_version == 1
	epsi=[];
	epsi(1,:)=1;                 % Zero Order (k=0) hermite
	epsi(2,:)=2*x;   % First Order (k=1) hermite
	if k==0
		y=epsi(1,:);
	elseif k==1
		y = epsi;
	else
		for n=3:k+1
			m=n-1;
			epsi(n,:)=(2*x*(epsi(n-1,:))-(2*(m-1)*(epsi(n-2,:))));
		end
		y = epsi;
    end
else
	epsi=[];
	epsi(1,:)=(pi)^(-1/4)*exp(-0.5*(x^2));                 % Zero Order (k=0) hermite
	epsi(2,:)=(((pi)^(-1/4)*exp(-0.5*(x^2)))*x)*sqrt(2);   % First Order (k=1) hermite
	if k==0
		y=epsi(1,:);
	elseif k==1
		y = epsi;
	else
		for n=3:k+1
			m=n-1;
			epsi(n,:)=(x*(sqrt(2/m))*(epsi(n-1,:)))-((sqrt((m-1)/m))*(epsi(n-2,:)));
		end
		y = epsi;
	end
	end
end