subroutine bdmk_eval(nd, ndim, eps, ikernel, beta, ipoly, norder, &
                        npbox, nboxes, nlevels, ltree, itree, iptr, &
                        centers, boxsize_in, fvals, pot)
  implicit none
  integer, intent(in) :: nd, ndim
  integer, intent(in) :: ikernel, ipoly, norder, npbox, nboxes, nlevels, ltree
  integer, intent(in) :: itree(ltree), iptr(8)
  real*8, intent(in) :: eps, beta
  real*8, intent(in) :: centers(ndim, nboxes)
  real*8, intent(in) :: boxsize_in(nlevels+1)
  real*8, intent(in) :: fvals(nd, npbox, nboxes)
  real*8, intent(inout) :: pot(nd, npbox, nboxes)

  integer :: jchnk, nd0, nhess, ifpgh, ifpghtarg, ntarg
  real*8, allocatable :: boxsize(:)
  real*8, allocatable :: fvalsj(:,:,:)
  real*8, allocatable :: pot0(:,:,:), grad0(:,:,:,:), hess0(:,:,:,:)
  real*8, allocatable :: targs(:,:), pote(:,:), grade(:,:,:), hesse(:,:,:)
  real*8, allocatable :: timeinfo(:)

  nd0 = 1
  nhess = ndim * (ndim + 1) / 2
  ifpgh = 1
  ifpghtarg = 0
  ntarg = 2

  allocate(boxsize(0:nlevels))
  boxsize = boxsize_in

  allocate(fvalsj(nd0, npbox, nboxes))
  allocate(pot0(nd0, npbox, nboxes))
  allocate(grad0(nd0, ndim, npbox, nboxes))
  allocate(hess0(nd0, nhess, npbox, nboxes))
  allocate(targs(ndim, ntarg))
  allocate(pote(nd0, ntarg))
  allocate(grade(nd0, ndim, ntarg))
  allocate(hesse(nd0, nhess, ntarg))
  allocate(timeinfo(20))

  do jchnk = 1, nd
    fvalsj = 0.0d0
    pot0 = 0.0d0
    grad0 = 0.0d0
    hess0 = 0.0d0
    targs = 0.0d0
    pote = 0.0d0
    grade = 0.0d0
    hesse = 0.0d0
    timeinfo = 0.0d0

    fvalsj(1,:,:) = fvals(jchnk,:,:)

    call bdmk(nd0, ndim, eps, ikernel, beta, ipoly, norder, npbox, &
              nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvalsj, &
              ifpgh, pot0, grad0, hess0, ntarg, targs, ifpghtarg, pote, &
              grade, hesse, timeinfo)

    pot(jchnk,:,:) = pot0(1,:,:)
  end do

  deallocate(timeinfo, hesse, grade, pote, targs)
  deallocate(hess0, grad0, pot0, fvalsj, boxsize)
end subroutine bdmk_eval
