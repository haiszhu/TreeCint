subroutine iddr_aid_mwrap(m, n, a, k, sk, rd, t)
  implicit none
  integer, intent(in) :: m, n, k
  double precision, intent(in) :: a(m, n)
  integer, intent(out) :: sk(k), rd(n-k)
  double precision, intent(out) :: t(k, n-k)

  integer :: i, j
  integer, allocatable :: list(:)
  double precision, allocatable :: proj(:), w(:)

  if (k < 1 .or. k >= n) then
    return
  end if

  allocate(list(n))
  allocate(proj(k*(n-k)))
  allocate(w((2*k+17)*n + 27*m + 100))

  call iddr_aidi(m, n, k, w)
  call iddr_aid(m, n, a, k, w, list, proj)

  do i = 1, k
    sk(i) = list(i)
  end do
  do i = 1, n-k
    rd(i) = list(k+i)
  end do

  do j = 1, n-k
    do i = 1, k
      t(i, j) = proj(i + k*(j-1))
    end do
  end do

  deallocate(w)
  deallocate(proj)
  deallocate(list)
end subroutine iddr_aid_mwrap
