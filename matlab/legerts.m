function [ts, whts] = legerts(itype,n,ts,whts)  % Gauss-Legendre nodes and weights on [-1,1]
% Tidied up from lgwt.m code by Greg von Winckel 02/25/2004, via Gillman.
% Included spectral diff matrix 10/23/2014, Bowei Wu
n=n-1; N1=n+1; N2=n+2;
xu=linspace(-1,1,N1)';
y=cos((2*(0:n)'+1)*pi/(2*n+2))+(0.27/N1)*sin(pi*xu*n/N2); % Initial guess
L=zeros(N1,N2);   % Legendre-Gauss Vandermonde Matrix
Lp=zeros(N1,N2);  % Derivative of LGVM
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2; % Iterate until new points are uniformly within epsilon of old points:
while max(abs(y-y0))>eps      
  L(:,1)=1; Lp(:,1)=0;
  L(:,2)=y; Lp(:,2)=1;
  for k=2:N1, L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k; end
  Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);       
  y0=y;
  y=y0-L(:,N2)./Lp;
end
ts = y(end:(-1):1);
whts=2./((1-y.^2).*Lp.^2)*(N2/N1)^2; % Compute the weights
end
