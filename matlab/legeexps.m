function [x, u, v, whts] = legeexps(itype, n, x, u, v, whts)

% Initialize the outputs
% x = zeros(n, 1);
% whts = zeros(n, 1);
% u = zeros(n, n);
% v = zeros(n, n);

% Determine the type of roots and weights to be calculated
itype_rts = 0;
if itype > 0
  itype_rts = 1;
end

% Call the legerts function to construct nodes and weights
% [x, whts] = legerts_mex(itype_rts, n, x, whts);
[x, whts] = legerts(itype_rts,n,x,whts);

% If itype is not 2, return early
if itype ~= 2
  return;
end

% Construct the matrix of values of the Legendre polynomials at these nodes
for i = 1:n
  % u(:, i) = legepols_mex(x(i), n-1, u(:, i));
  u(:, i) = legepols(x(i), n-1, u(:, i));
end

% Construct the matrix v by transposing u
v = u.';

% Construct the inverse u by converting values at Gaussian nodes into coefficients
for i = 1:n
  d = 1;
  d = d * (2 * i - 1) / 2;
  for j = 1:n
    u(i, j) = v(j, i) * whts(j) * d;
  end
end
end