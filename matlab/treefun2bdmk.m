function [src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly)

dpars = zeros(1000,1);
ipars = zeros(100,1);
zpars = zeros(10,1);
% dada
dom = ratio*f.domain(:,1)';
scale = (dom(2)-dom(1))/2;
norder = f.n;
nhess = ndim*(ndim+1)/2;
% convert to dmk tree format (should not change the tree if opts.balance == true)
iper = 0;
npbox = norder^ndim;
grid = zeros(ndim,norder^2);
nbmax = 2^ndim*f.id(end); % is 2^ndim enough, since some box needs to be refined more than once...
nlmax = f.height(1); % is this ok?
centers0 = ratio/2*(f.domain(1:2:end,:)+f.domain(2:2:end,:));
centers = zeros(ndim,nbmax); centers(:,1:size(centers0,2)) = centers0;
nlevels = f.height(1);
nboxes = f.id(end);
boxsize = 2*scale./2.^(0:nlmax); 
laddr = zeros(2,nlmax+1); 
for ilev=0:nlmax
  laddr(1,ilev+1) = sum(f.level<ilev)+1;
  laddr(2,ilev+1) = sum(f.level<=ilev);
end
ilevel = zeros(nbmax,1); ilevel(1:nboxes) = f.level;
iparent = zeros(nbmax,1); iparent(1:nboxes) = f.parent;
nchild = zeros(nbmax,1); nchild(1:nboxes) = sum(f.children>0);
ichild = zeros(2^ndim,nbmax); ichild(:,1:nboxes) = f.children; ichild(ichild==0) = -1; % -1 convention in fortran tree
% this does not belong to fortran tree, but needed for treefun plot
coeffs = cell(1,nbmax); coeffs(1:nboxes) = f.coeffs(1:nboxes);
% call computecoll to get vol_tree_fix_lr coll info
nnbors0 = zeros(nboxes,1); nbors0 = zeros(3^ndim,nboxes);
[nnbors0,nbors0] = ...
  computecoll(ndim,nlevels,nboxes,laddr,boxsize,...
                  centers(:,1:nboxes),iparent(1:nboxes),nchild(1:nboxes),ichild(:,1:nboxes),iper,...
                  nnbors0,nbors0);
nnbors = zeros(nbmax,1); nbors = zeros(3^ndim,nbmax);
nnbors(1:nboxes) = nnbors0; nbors(:,1:nboxes) = nbors0;
% now fix lr (with coeffs, which is different from)
ndtmp = 1;
[centers,nboxes,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,coeffs] = ...
      vol_tree_fix_lr_treefun(ndim,iper,ndtmp,dpars,zpars,ipars,...
                               norder,npbox,grid,centers,nlevels,nboxes,boxsize,...
                               nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,...
                               coeffs);
centers = centers(:,1:nboxes);

% after fix lr
iperiod=0;
itype=0;
zk=60.0d0;
mc=2^ndim;
mnbors=3^ndim;
npbox=norder^ndim;
nleafbox=numel(f.leaves);
boxlen=boxsize(1);
iptype=2;
eta=0.0d0;
rintl=zeros(nlevels+1,1);
iptr=zeros(8,1);
iptr(1) = 1;
iptr(2) = 2*(nlevels+1)+1;
iptr(3) = iptr(2) + nboxes;
iptr(4) = iptr(3) + nboxes;
iptr(5) = iptr(4) + nboxes;
iptr(6) = iptr(5) + mc*nboxes;
iptr(7) = iptr(6) + nboxes;
iptr(8) = iptr(7) + mnbors*nboxes;
ltree = (4 + mc + mnbors) * nboxes + 2 * (nlevels + 1);
%
itree = zeros(ltree,1);
itree(iptr(1):(iptr(1)+2*nlevels+1)) = reshape(laddr,[],1);
itree(iptr(2):(iptr(2)+nboxes-1)) = ilevel(1:nboxes);
itree(iptr(3):(iptr(3)+nboxes-1)) = iparent(1:nboxes);
itree(iptr(4):(iptr(4)+nboxes-1)) = nchild(1:nboxes);
itree(iptr(5):(iptr(5)+2^ndim*nboxes-1)) = reshape(ichild(:,1:nboxes),[],1);
itree(iptr(6):(iptr(6)+nboxes-1)) = nnbors(1:nboxes);
itree(iptr(7):(iptr(7)+3^ndim*nboxes-1)) = reshape(nbors(:,1:nboxes),[],1);
%
src = zeros(ndim,npbox,nboxes);
wts = zeros(npbox,nboxes);
srcleaf = zeros(ndim,npbox,nleafbox);
wtsleaf = zeros(npbox,nleafbox);
xq = zeros(norder, 1);
wtsq = zeros(norder, 1);
umat = zeros(norder, norder);
vmat = zeros(norder, norder);
grid = zeros(ndim, npbox);
wtsb = zeros(npbox, 1);
if ipoly == 0
  [xq, umat, vmat, wtsq] = legeexps(itype, norder, xq, umat, vmat, wtsq);
elseif ipoly == 1
  [xq, umat, vmat, wtsq] = chebexps(itype, norder, xq, umat, vmat, wtsq);
end
xq = xq / 2;
wtsq = wtsq / 2;
grid = mesh3d(xq, norder, xq, norder, xq, norder, grid);
ind = 0;
for iz = 1:norder
  for iy = 1:norder
    for ix = 1:norder
      ind = ind + 1;
      wtsb(ind) = wtsq(ix) * wtsq(iy) * wtsq(iz);
    end
  end
end
% Main loop to compute src on the grid points
jbox = 0;
xyz = zeros(ndim,1);
for ilev = 0:nlevels
  bs = boxsize(ilev+1); % Adjusted for MATLAB indexing (1-based)
  ifirstbox = itree(2 * ilev + 1);
  ilastbox = itree(2 * ilev + 2);
  nbloc = ilastbox - ifirstbox + 1;
  for i = 1:nbloc
    ibox = ifirstbox + i - 1;
    for ell = 1:npbox
      for j = 1:ndim
        xyz(j) = centers(j, ibox) + grid(j, ell) * bs;
        src(j, ell, ibox) = xyz(j);
      end
      wts(ell, ibox) = wtsb(ell) * bs^ndim;
    end
    if itree(iptr(4) + ibox - 1) == 0
      jbox = jbox + 1;
      for ell = 1:npbox
        for j = 1:ndim
          srcleaf(j, ell, jbox) = src(j, ell, ibox);
        end
        wtsleaf(ell, jbox) = wts(ell, ibox);
      end
    end
  end
end

end