function [SK, RD, T] = iddr_aid_mwrap_mex(A, krank, SK, RD, T)
[m, n] = size(A);
k = min([floor(krank), m, n]);
k = max(k, 1);
if k >= n
  error('iddr_aid_mwrap_mex requires 1 <= krank < size(A,2).');
end
nrd = n - k;
mex_id_ = 'iddr_aid_mwrap(c i int[x], c i int[x], c i double[xx], c i int[x], c io int[x], c io int[x], c io double[xx])';
[SK, RD, T] = libid_mex(mex_id_, m, n, A, k, SK, RD, T, 1, 1, m, n, 1, k, nrd, k, nrd);
end
