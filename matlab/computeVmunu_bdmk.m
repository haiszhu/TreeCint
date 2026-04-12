function Vmunu = computeVmunu_bdmk(fvals, nleafbox, wtsleaf, ratio, ...
    ndim, eps, ikernel, beta, ipoly, norder, npbox, nboxes, nlevels, ...
    ltree, itree, iptr, centers, boxsize)
% Compute Vmunu from interpolating vectors using in-process BDMK MEX.
% fvals is (npts x nd) where npts = npbox*nboxes.

nd = size(fvals,2);
phi_kl = zeros(nd,npbox,nboxes);
for k = 1:nd
  phi_kl(k,:,:) = reshape(fvals(:,k),[npbox nboxes]);
end

pot = zeros(nd,npbox,nboxes);
pot = bdmk_eval_mex(nd, ndim, eps, ikernel, beta, ipoly, norder, npbox, ...
  nboxes, nlevels, ltree, itree, iptr, centers, boxsize, phi_kl, pot);
pot = reshape(pot, [nd, npbox, nboxes]);

phi_ij_leaf = zeros(nd,npbox,nleafbox);
potleaf = zeros(nd,npbox,nleafbox);
jbox = 0;
for ilev = 0:nlevels
  ifirstbox = itree(2 * ilev + 1);
  ilastbox = itree(2 * ilev + 2);
  nbloc = ilastbox - ifirstbox + 1;
  for i = 1:nbloc
    ibox = ifirstbox + i - 1;
    if itree(iptr(4) + ibox - 1) == 0
      jbox = jbox + 1;
      phi_ij_leaf(:,:,jbox) = phi_kl(:,:,ibox);
      potleaf(:,:,jbox) = pot(:,:,ibox);
    end
  end
end

potleaf = reshape(potleaf, nd, []) / ratio^2;
phi_ij_leaf = reshape(phi_ij_leaf, nd, []);
wtsleaf = wtsleaf(:)';
phi_ij_leaf = phi_ij_leaf .* wtsleaf;
Vmunu = (phi_ij_leaf * potleaf') / ratio^3;
end
