! bdmk_scf.f90 — Standalone Fortran program: ISDF → Vmunu via Box-DMK
!
! Usage:
!   ./int2-bdmk-scf <treefun.h5> <isdf.h5> <output.h5>
!
! Reads tree metadata and ISDF interpolating vectors from HDF5, runs the
! Box-DMK Coulomb solver (OMP-parallel over nd components), assembles Vmunu
! via a single DGEMM, and writes Vmunu to HDF5.
!
! Build (see Makefile target int2-bdmk-scf):
!   gfortran -O3 -fopenmp -fallow-argument-mismatch -w \
!     -I/opt/homebrew/opt/hdf5/include \
!     int2_bdmk_scf.f90 libbdmk_ref.a \
!     -L/opt/homebrew/opt/hdf5/lib -lhdf5_fortran -lhdf5 \
!     -framework Accelerate -lm -o int2-bdmk-scf
!
! Design notes:
!   - interp_vecs(1,:,k) is read contiguously inside the OMP loop (no fvals
!     staging array needed); avoids the strided triple-loop of the reference.
!   - Scratch arrays (fvalsj, pot0, ...) allocated once per thread outside
!     the !$omp do, not once per iteration.
!   - Vmunu assembled with a single dgemm('N','T',...) call (level-3 BLAS)
!     rather than nd separate MATVEC calls.
!   - beta_k holds the kernel parameter; alpha_d/beta_d are DGEMM scalars
!     (no variable shadowing).
!   - dims2 is passed consistently to both h5screate_simple_f and h5dwrite_f
!     for the Vmunu dataset.

      program int2_bdmk_scf
      use omp_lib
      use hdf5
      implicit none

      ! ----- HDF5 handles -----
      integer(HID_T) :: fid_tree, fid_isdf, fid_out
      integer(HID_T) :: dset_id, dspace_id
      integer(HSIZE_T) :: dims1(1), dims2(2)
      integer(HSIZE_T) :: isdf_dims(3), isdf_maxdims(3)
      integer :: herr

      ! ----- Tree metadata -----
      integer :: ndim, ikernel, ipoly, norder, npbox, nboxes
      integer :: nlevels, ltree, nleafbox
      real*8  :: eps, beta_k, ratio       ! beta_k: kernel param (not DGEMM scalar)
      integer, allocatable :: itree(:)
      integer :: iptr(8)
      real*8,  allocatable :: centers(:,:)   ! (ndim, nboxes)
      real*8,  allocatable :: boxsize(:)     ! (0:nlevels)
      real*8,  allocatable :: wtsleaf(:,:)   ! (npbox, nleafbox)

      ! ----- ISDF data -----
      integer :: nd, npts
      real*8, allocatable :: interp_vecs(:,:,:)  ! (2, npts, nd)

      ! ----- BDMK output -----
      real*8, allocatable :: pot(:,:,:)   ! (nd, npbox, nboxes)

      ! ----- OMP private scratch (declared here; allocated inside parallel) -----
      real*8, allocatable :: fvalsj(:,:,:)      ! (1, npbox, nboxes)
      real*8, allocatable :: pot0(:,:,:)         ! (1, npbox, nboxes)
      real*8, allocatable :: grad0(:,:,:,:)
      real*8, allocatable :: hess0(:,:,:,:)
      real*8, allocatable :: targs(:,:)
      real*8, allocatable :: pote(:,:)
      real*8, allocatable :: grade(:,:,:)
      real*8, allocatable :: hesse(:,:,:)
      real*8, allocatable :: timeinfo(:)

      ! ----- Leaf extraction -----
      integer :: nleafbox_actual
      integer, allocatable :: iboxleaf(:)
      real*8, allocatable :: phi_leaf(:,:)   ! (nd, npbox*nleafbox), Ask at leaves
      real*8, allocatable :: pot_leaf(:,:)   ! (nd, npbox*nleafbox), pot at leaves

      ! ----- Vmunu -----
      real*8, allocatable :: Vmunu(:,:)   ! (nd, nd)
      real*8 :: alpha_d, beta_d           ! DGEMM scalars (distinct from beta_k)

      ! ----- Loop/misc -----
      integer :: i, k, jchnk, jl, ibox
      integer :: ilev, ifirstbox, ilastbox, nbloc
      integer :: nd0, nhess, ifpgh, ifpghtarg, ntarg
      real*8  :: t1, t2, t_bdmk, t_vmunu, t_total

      character(len=256) :: treefun_filename, isdf_filename, output_filename
      integer :: nargs

      ! Fixed BDMK call parameters
      nd0        = 1
      ifpgh      = 1
      ifpghtarg  = 0
      ntarg      = 2

      ! ----- CLI arguments -----
      nargs = command_argument_count()
      if (nargs >= 3) then
        call get_command_argument(1, treefun_filename)
        call get_command_argument(2, isdf_filename)
        call get_command_argument(3, output_filename)
      else
        treefun_filename = 'treefun_h2o_cc-pvdz.h5'
        isdf_filename    = 'myisdf_1e-3.h5'
        output_filename  = 'bdmk_1e-3.h5'
      endif
      print *, 'treefun : ', trim(treefun_filename)
      print *, 'isdf    : ', trim(isdf_filename)
      print *, 'output  : ', trim(output_filename)
      t_total = omp_get_wtime()

      ! ===== Initialize HDF5 =====
      call h5open_f(herr)

      ! ===== Read treefun HDF5 =====
      call h5fopen_f(treefun_filename, H5F_ACC_RDONLY_F, fid_tree, herr)
      if (herr /= 0) then; print *, 'ERROR: cannot open ', trim(treefun_filename); stop; endif

      dims1 = [1]
      call h5dopen_f(fid_tree, '/ndim',     dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ndim,     dims1, herr)
      call h5dopen_f(fid_tree, '/ikernel',  dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ikernel,  dims1, herr)
      call h5dopen_f(fid_tree, '/ipoly',    dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ipoly,    dims1, herr)
      call h5dopen_f(fid_tree, '/norder',   dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, norder,   dims1, herr)
      call h5dopen_f(fid_tree, '/npbox',    dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, npbox,    dims1, herr)
      call h5dopen_f(fid_tree, '/nboxes',   dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nboxes,   dims1, herr)
      call h5dopen_f(fid_tree, '/nlevels',  dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nlevels,  dims1, herr)
      call h5dopen_f(fid_tree, '/ltree',    dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ltree,    dims1, herr)
      call h5dopen_f(fid_tree, '/nleafbox', dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nleafbox, dims1, herr)
      call h5dopen_f(fid_tree, '/eps',      dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, eps,       dims1, herr)
      call h5dopen_f(fid_tree, '/beta',     dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, beta_k,    dims1, herr)
      call h5dopen_f(fid_tree, '/ratio',    dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ratio,     dims1, herr)

      allocate(itree(ltree))
      dims1 = [ltree]
      call h5dopen_f(fid_tree, '/itree', dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, itree, dims1, herr)

      dims1 = [8]
      call h5dopen_f(fid_tree, '/iptr', dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_INTEGER, iptr, dims1, herr)

      allocate(centers(ndim, nboxes))
      dims2 = [ndim, nboxes]
      call h5dopen_f(fid_tree, '/centers', dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, centers, dims2, herr)

      allocate(boxsize(0:nlevels))
      dims1 = [nlevels+1]
      call h5dopen_f(fid_tree, '/boxsize', dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, boxsize, dims1, herr)

      allocate(wtsleaf(npbox, nleafbox))
      dims2 = [npbox, nleafbox]
      call h5dopen_f(fid_tree, '/wtsleaf', dset_id, herr)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wtsleaf, dims2, herr)

      call h5dclose_f(dset_id, herr)
      call h5fclose_f(fid_tree, herr)

      print *, 'ndim=', ndim, ' norder=', norder, ' npbox=', npbox, &
               ' nboxes=', nboxes, ' nlevels=', nlevels, ' nleafbox=', nleafbox
      print *, 'eps=', eps, ' beta=', beta_k, ' ikernel=', ikernel, ' ratio=', ratio

      ! ===== Read ISDF HDF5 =====
      call h5fopen_f(isdf_filename, H5F_ACC_RDONLY_F, fid_isdf, herr)
      if (herr /= 0) then; print *, 'ERROR: cannot open ', trim(isdf_filename); stop; endif

      call h5dopen_f(fid_isdf, '/interpolating_vectors', dset_id, herr)
      call h5dget_space_f(dset_id, dspace_id, herr)
      call h5sget_simple_extent_dims_f(dspace_id, isdf_dims, isdf_maxdims, herr)
      ! isdf_dims = (2, npts, nd)  [Fortran column-major matches MATLAB write]
      nd   = isdf_dims(3)
      npts = isdf_dims(2)
      allocate(interp_vecs(isdf_dims(1), isdf_dims(2), isdf_dims(3)))
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, interp_vecs, isdf_dims, herr)
      call h5dclose_f(dset_id, herr)
      call h5sclose_f(dspace_id, herr)
      call h5fclose_f(fid_isdf, herr)

      print *, 'nd=', nd, ' npts=', npts
      if (npts /= npbox*nboxes) then
        print *, 'ERROR: npts mismatch: expected', npbox*nboxes, 'got', npts; stop
      endif

      ! ===== OMP-parallel BDMK =====
      ! Key design: read interp_vecs(1,:,jchnk) directly in each OMP iteration.
      ! interp_vecs(1,:,k) is contiguous in Fortran column-major for fixed k,
      ! avoiding the strided fvals(k,:,:) gather of the naive approach.
      !
      ! Scratch arrays allocated once per thread (outside !$omp do) so the
      ! allocator is called nthreads times, not nd times.

      nhess = ndim * (ndim+1) / 2
      allocate(pot(nd, npbox, nboxes))
      pot = 0.0d0

      print *, '=========Start BDMK======='
      t1 = omp_get_wtime()

      !$omp parallel default(none) &
      !$omp& shared(nd, ndim, eps, ikernel, beta_k, ipoly, norder, npbox, nboxes, &
      !$omp&        nlevels, ltree, itree, iptr, centers, boxsize, &
      !$omp&        interp_vecs, pot, nd0, nhess, ifpgh, ifpghtarg, ntarg) &
      !$omp& private(jchnk, fvalsj, pot0, grad0, hess0, targs, pote, grade, hesse, timeinfo)

        allocate(fvalsj  (nd0, npbox, nboxes))
        allocate(pot0    (nd0, npbox, nboxes))
        allocate(grad0   (nd0, ndim,  npbox, nboxes))
        allocate(hess0   (nd0, nhess, npbox, nboxes))
        allocate(targs   (ndim, ntarg))
        allocate(pote    (nd0, ntarg))
        allocate(grade   (nd0, ndim,  ntarg))
        allocate(hesse   (nd0, nhess, ntarg))
        allocate(timeinfo(20))

        !$omp do schedule(static)
        do jchnk = 1, nd
          fvalsj   = 0.0d0
          pot0     = 0.0d0
          grad0    = 0.0d0
          hess0    = 0.0d0
          targs    = 0.0d0
          pote     = 0.0d0
          grade    = 0.0d0
          hesse    = 0.0d0
          timeinfo = 0.0d0

          ! Contiguous read: interp_vecs(1,:,jchnk) is a npts-element
          ! column-major slice (first index fixed, second varies fastest)
          fvalsj(1,:,:) = reshape(interp_vecs(1,:,jchnk), [npbox, nboxes])

          call bdmk(nd0, ndim, eps, ikernel, beta_k, ipoly, norder, npbox, &
                    nboxes, nlevels, ltree, itree, iptr, centers, boxsize,  &
                    fvalsj, ifpgh, pot0, grad0, hess0, ntarg, targs,        &
                    ifpghtarg, pote, grade, hesse, timeinfo)

          pot(jchnk,:,:) = pot0(1,:,:)
        enddo
        !$omp end do

        deallocate(fvalsj, pot0, grad0, hess0)
        deallocate(targs, pote, grade, hesse, timeinfo)

      !$omp end parallel

      t2 = omp_get_wtime()
      t_bdmk = t2 - t1
      print *, '=========End BDMK======='
      print *, 'BDMK time:', t_bdmk, 'seconds with', omp_get_max_threads(), 'threads'

      ! ===== Precompute leaf box indices =====
      allocate(iboxleaf(nleafbox))
      nleafbox_actual = 0
      do ilev = 0, nlevels
        ifirstbox = itree(2*ilev + 1)
        ilastbox  = itree(2*ilev + 2)
        nbloc = ilastbox - ifirstbox + 1
        do i = 1, nbloc
          ibox = ifirstbox + i - 1
          if (itree(iptr(4) + ibox - 1) == 0) then
            nleafbox_actual = nleafbox_actual + 1
            iboxleaf(nleafbox_actual) = ibox
          endif
        enddo
      enddo
      if (nleafbox_actual /= nleafbox) then
        print *, 'WARNING: leaf count mismatch:', nleafbox_actual, 'vs', nleafbox
      endif

      ! ===== Extract Ask and pot at leaf boxes =====
      ! phi_leaf(nd, npbox*nleafbox): Ask values at leaf quadrature points
      ! pot_leaf(nd, npbox*nleafbox): BDMK potential, scaled by 1/ratio^2
      !
      ! pot(:,:,ibox) is a contiguous (nd*npbox) block in Fortran column-major.
      ! interp_vecs(1,:,k) accessed as a contiguous slice for fixed k.
      allocate(phi_leaf(nd, npbox * nleafbox))
      allocate(pot_leaf(nd, npbox * nleafbox))

      do jl = 1, nleafbox
        ibox = iboxleaf(jl)
        ! pot(:,:,ibox) is contiguous: nd*npbox elements
        pot_leaf(:, (jl-1)*npbox+1 : jl*npbox) = pot(:, :, ibox) / ratio**2
        ! Build phi_leaf from interp_vecs row by row (each row k is contiguous)
        do k = 1, nd
          phi_leaf(k, (jl-1)*npbox+1 : jl*npbox) = &
            interp_vecs(1, (ibox-1)*npbox+1 : ibox*npbox, k)
        enddo
      enddo

      deallocate(pot, iboxleaf, interp_vecs)

      ! Apply quadrature weights to phi_leaf in-place: phi_leaf(:,r) *= wts(r)
      do jl = 1, nleafbox
        do i = 1, npbox
          phi_leaf(:, (jl-1)*npbox + i) = &
            phi_leaf(:, (jl-1)*npbox + i) * wtsleaf(i, jl)
        enddo
      enddo
      deallocate(wtsleaf)

      ! ===== Vmunu via single DGEMM =====
      ! Vmunu(mu,nu) = sum_r phi_leaf_weighted(mu,r) * pot_leaf(nu,r) / ratio^3
      !              = (phi_leaf * pot_leaf^T)(mu,nu) / ratio^3
      !
      ! dgemm('N','T', nd, nd, nleaf_pts, alpha, phi_leaf, nd, pot_leaf, nd, beta, Vmunu, nd)
      ! where alpha = 1/ratio^3, beta = 0.
      !
      ! This is a single level-3 BLAS call (compute-bound, cache-blocked by BLAS),
      ! replacing nd separate MATVEC calls.
      allocate(Vmunu(nd, nd))

      print *, '=========Start Vmunu======='
      t1 = omp_get_wtime()

      alpha_d = 1.0d0 / ratio**3
      beta_d  = 0.0d0
      call dgemm('N', 'T', nd, nd, npbox*nleafbox, &
                 alpha_d, phi_leaf, nd, pot_leaf, nd, &
                 beta_d,  Vmunu,    nd)

      t2 = omp_get_wtime()
      t_vmunu = t2 - t1
      print *, '=========End Vmunu======='
      print *, 'Vmunu time:', t_vmunu, 'seconds'

      deallocate(phi_leaf, pot_leaf)

      ! ===== Write Vmunu to output HDF5 =====
      call h5fcreate_f(output_filename, H5F_ACC_TRUNC_F, fid_out, herr)
      if (herr /= 0) then; print *, 'ERROR: cannot create ', trim(output_filename); stop; endif

      dims2 = [nd, nd]
      call h5screate_simple_f(2, dims2, dspace_id, herr)
      call h5dcreate_f(fid_out, '/Vmunu', H5T_NATIVE_DOUBLE, dspace_id, dset_id, herr)
      call h5dwrite_f (dset_id, H5T_NATIVE_DOUBLE, Vmunu, dims2, herr)
      call h5dclose_f(dset_id, herr)
      call h5sclose_f(dspace_id, herr)
      call h5fclose_f(fid_out, herr)

      call h5close_f(herr)

      t_total = omp_get_wtime() - t_total
      print *, 'Wrote Vmunu [', nd, 'x', nd, '] to: ', trim(output_filename)
      print *, 'Total time:', t_total, 'seconds'

      deallocate(Vmunu, itree, centers, boxsize)

      end program int2_bdmk_scf
