program main
  use baobzi
  implicit none
  type(c_funptr) :: func
  real(kind=c_double) :: center(2), half_length(2), tol
  real(kind=c_double) :: x(2)
  integer(kind=c_int16_t) :: order, dim
  type(c_ptr) :: func_approx
  character(len=64) :: fname

  func = c_funloc(testfun)
  center(:) = 0.0
  half_length(:) = 1.0
  tol = 1E-8
  dim = 2
  order = 6

  x(:) = 0.25

  fname = trim(adjustl('fortran.baobzi'))//char(0)

  func_approx = baobzi_init(func, dim, order, center, half_length, tol)
  print *, baobzi_eval(func_approx, x) - testfun(x)

  call baobzi_save(func_approx, fname)
  func_approx = baobzi_free(func_approx)

  func_approx = baobzi_restore(func, fname)
  print *, baobzi_eval(func_approx, x) - testfun(x)

contains
  function testfun (x) bind(c) result(y)
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=c_double), dimension(*) :: x
    real(kind=c_double) :: y

    y = x(1) * x(2)
  end function testfun

end program main
