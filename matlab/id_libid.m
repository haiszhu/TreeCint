function [SK, RD, T] = id_libid(A, k)
%ID_LIBID Rank-based interpolative decomposition via libid (iddr_aid).
%   [SK,RD,T] = ID_LIBID(A,K) returns skeleton indices SK, redundant
%   indices RD, and interpolation matrix T such that
%   A(:,RD) ~= A(:,SK) * T.

[m, n] = size(A);
kk = min([floor(k), m, n]);
if kk < 1 || kk >= n
  error('id_libid:invalidRank', 'Require 1 <= k < size(A,2).');
end
[SK, RD, T] = iddr_aid_mwrap_mex(A, kk, zeros(1, kk), zeros(1, n-kk), zeros(kk, n-kk));
end
