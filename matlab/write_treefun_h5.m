function write_treefun_h5(molname,basmod,ndim,treefun_eps,ikernel,beta,...
                          ipoly,norder,npbox,nboxes,nlevels,ltree,itree,...
                          iptr,centers,boxsize,DS1,nleafbox,wtsleaf,ratio)
%
%

[~,bas_string,~] = fileparts(basmod); 
treefun_h5_fn = ['treefun_' molname '_' bas_string '.h5']; % treefun .h5 filename
file_id = H5F.create(treefun_h5_fn, "H5F_ACC_TRUNC", "H5P_DEFAULT", "H5P_DEFAULT"); % Overwrites any existing file with the same name
% save ndim
ndim_dataspace_id = H5S.create_simple(1, 1, []);
ndim_dcpl_id = H5P.create('H5P_DATASET_CREATE');
ndim_dataset_id = H5D.create(file_id, '/ndim', 'H5T_NATIVE_INT', ndim_dataspace_id, ndim_dcpl_id);
H5D.write(ndim_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(ndim));
% save treefun_eps
eps_dataspace_id = H5S.create_simple(1, 1, []);
eps_dcpl_id = H5P.create('H5P_DATASET_CREATE');
eps_dataset_id = H5D.create(file_id, '/eps', 'H5T_NATIVE_DOUBLE', eps_dataspace_id, eps_dcpl_id);
H5D.write(eps_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', treefun_eps);
% save ikernel
ikernel_dataspace_id = H5S.create_simple(1, 1, []);
ikernel_dcpl_id = H5P.create('H5P_DATASET_CREATE');
ikernel_dataset_id = H5D.create(file_id, '/ikernel', 'H5T_NATIVE_INT', ikernel_dataspace_id, ikernel_dcpl_id);
H5D.write(ikernel_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(ikernel));
% save beta
beta_dataspace_id = H5S.create_simple(1, 1, []);
beta_dcpl_id = H5P.create('H5P_DATASET_CREATE');
beta_dataset_id = H5D.create(file_id, '/beta', 'H5T_NATIVE_DOUBLE', beta_dataspace_id, beta_dcpl_id);
H5D.write(beta_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', beta);
% save ipoly
ipoly_dataspace_id = H5S.create_simple(1, 1, []);
ipoly_dcpl_id = H5P.create('H5P_DATASET_CREATE');
ipoly_dataset_id = H5D.create(file_id, '/ipoly', 'H5T_NATIVE_INT', ipoly_dataspace_id, ipoly_dcpl_id);
H5D.write(ipoly_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(ipoly));
% save norder
norder_dataspace_id = H5S.create_simple(1, 1, []);
norder_dcpl_id = H5P.create('H5P_DATASET_CREATE');
norder_dataset_id = H5D.create(file_id, '/norder', 'H5T_NATIVE_INT', norder_dataspace_id, norder_dcpl_id);
H5D.write(norder_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(norder));
% save npbox
npbox_dataspace_id = H5S.create_simple(1, 1, []);
npbox_dcpl_id = H5P.create('H5P_DATASET_CREATE');
npbox_dataset_id = H5D.create(file_id, '/npbox', 'H5T_NATIVE_INT', npbox_dataspace_id, npbox_dcpl_id);
H5D.write(npbox_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(npbox));
% save nboxes
nboxes_dataspace_id = H5S.create_simple(1, 1, []);
nboxes_dcpl_id = H5P.create('H5P_DATASET_CREATE');
nboxes_dataset_id = H5D.create(file_id, '/nboxes', 'H5T_NATIVE_INT', nboxes_dataspace_id, nboxes_dcpl_id);
H5D.write(nboxes_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(nboxes));
% save nlevels
nlevels_dataspace_id = H5S.create_simple(1, 1, []);
nlevels_dcpl_id = H5P.create('H5P_DATASET_CREATE');
nlevels_dataset_id = H5D.create(file_id, '/nlevels', 'H5T_NATIVE_INT', nlevels_dataspace_id, nlevels_dcpl_id);
H5D.write(nlevels_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(nlevels));
% save ltree
ltree_dataspace_id = H5S.create_simple(1, 1, []);
ltree_dcpl_id = H5P.create('H5P_DATASET_CREATE');
ltree_dataset_id = H5D.create(file_id, '/ltree', 'H5T_NATIVE_INT', ltree_dataspace_id, ltree_dcpl_id);
H5D.write(ltree_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(ltree));
% save itree
itree_dataspace_id = H5S.create_simple(1, ltree, []);
itree_dcpl_id = H5P.create('H5P_DATASET_CREATE');
itree_dataset_id = H5D.create(file_id, '/itree', 'H5T_NATIVE_INT', itree_dataspace_id, itree_dcpl_id);
H5D.write(itree_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(itree));
% save iptr
iptr_dataspace_id = H5S.create_simple(1, 8, []);
iptr_dcpl_id = H5P.create('H5P_DATASET_CREATE');
iptr_dataset_id = H5D.create(file_id, '/iptr', 'H5T_NATIVE_INT', iptr_dataspace_id, iptr_dcpl_id);
H5D.write(iptr_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(iptr));
% save centers
centers_dataspace_id = H5S.create_simple(2, fliplr(size(centers)), []);
centers_dcpl_id = H5P.create('H5P_DATASET_CREATE');
centers_dataset_id = H5D.create(file_id, '/centers', 'H5T_NATIVE_DOUBLE', centers_dataspace_id, centers_dcpl_id);
H5D.write(centers_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', centers);
% save boxsize
boxsize_dataspace_id = H5S.create_simple(1, nlevels+1, []);
boxsize_dcpl_id = H5P.create('H5P_DATASET_CREATE');
boxsize_dataset_id = H5D.create(file_id, '/boxsize', 'H5T_NATIVE_DOUBLE', boxsize_dataspace_id, boxsize_dcpl_id);
H5D.write(boxsize_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', boxsize);
% save src
src_dataspace_id = H5S.create_simple(3, fliplr(size(DS1)), []);
src_dcpl_id = H5P.create('H5P_DATASET_CREATE');
src_dataset_id = H5D.create(file_id, '/DS1', 'H5T_NATIVE_DOUBLE', src_dataspace_id, src_dcpl_id);
H5D.write(src_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', DS1);
% save nleafbox
nleafbox_dataspace_id = H5S.create_simple(1, 1, []);
nleafbox_dcpl_id = H5P.create('H5P_DATASET_CREATE');
nleafbox_dataset_id = H5D.create(file_id, '/nleafbox', 'H5T_NATIVE_INT', nleafbox_dataspace_id, nleafbox_dcpl_id);
H5D.write(nleafbox_dataset_id, 'H5T_NATIVE_INT32', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', int32(nleafbox));
% save wtsleaf
wtsleaf_dataspace_id = H5S.create_simple(2, fliplr(size(wtsleaf)), []);
wtsleaf_dcpl_id = H5P.create('H5P_DATASET_CREATE');
wtsleaf_dataset_id = H5D.create(file_id, '/wtsleaf', 'H5T_NATIVE_DOUBLE', wtsleaf_dataspace_id, wtsleaf_dcpl_id);
H5D.write(wtsleaf_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', wtsleaf);
% save ratio
ratio_dataspace_id = H5S.create_simple(1, 1, []);
ratio_dcpl_id = H5P.create('H5P_DATASET_CREATE');
ratio_dataset_id = H5D.create(file_id, '/ratio', 'H5T_NATIVE_DOUBLE', ratio_dataspace_id, ratio_dcpl_id);
H5D.write(ratio_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', ratio);


H5D.close(ndim_dataset_id);
H5D.close(eps_dataset_id);
H5D.close(ikernel_dataset_id);
H5D.close(beta_dataset_id);
H5D.close(ipoly_dataset_id);
H5D.close(norder_dataset_id);
H5D.close(npbox_dataset_id);
H5D.close(nboxes_dataset_id);
H5D.close(nlevels_dataset_id);
H5D.close(ltree_dataset_id);
H5D.close(itree_dataset_id);
H5D.close(iptr_dataset_id);
H5D.close(centers_dataset_id);
H5D.close(boxsize_dataset_id);
H5D.close(src_dataset_id);
H5D.close(nleafbox_dataset_id);
H5D.close(wtsleaf_dataset_id);
H5D.close(ratio_dataset_id);
H5F.close(file_id);

end