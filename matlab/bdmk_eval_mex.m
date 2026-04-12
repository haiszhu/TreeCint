function [pot] = bdmk_eval_mex(nd, ndim, eps, ikernel, beta, ipoly, norder, npbox, nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals, pot)
nboxsize = numel(boxsize);
nfvals = numel(fvals);
npot = numel(pot);
mex_id_ = 'bdmk_eval(c i int[x], c i int[x], c i double[x], c i int[x], c i double[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i int[x], c i double[xx], c i double[x], c i double[x], c io double[x])';
[pot] = bdmk_mex(mex_id_, nd, ndim, eps, ikernel, beta, ipoly, norder, npbox, nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals, pot, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ltree, 8, ndim, nboxes, nboxsize, nfvals, npot);
end
