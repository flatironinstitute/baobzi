module baobzi
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function baobzi_init (fin, dim, order, center, half_length, tol) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      type(c_funptr), intent(in), value :: fin
      integer(kind=c_int16_t), intent(in), value :: dim, order
      real(kind=c_double), dimension(*) :: center, half_length
      real(kind=c_double), intent(in), value :: tol
      type(c_ptr) :: func
    end function baobzi_init

    subroutine baobzi_save (fin, fname) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: fin
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_ptr) :: func
    end subroutine baobzi_save

    function baobzi_restore (f, fname) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      type(c_funptr), intent(in), value :: f
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_ptr) :: func
    end function baobzi_restore

    function baobzi_eval (f, x) bind(c) result(y)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: f
      real(kind=c_double), dimension(*) :: x
      real(kind=c_double) :: y
    end function baobzi_eval

    function baobzi_free (f) bind(c) result(f_null)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: f
      type(c_ptr) :: f_null
    end function baobzi_free
  end interface

end module baobzi
