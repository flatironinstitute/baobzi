module baobzi
  use, intrinsic :: iso_c_binding
  implicit none

  type, bind(c) :: baobzi_input_t
    type(c_funptr) :: func = c_null_funptr
    type(c_ptr) :: data = c_null_ptr
    integer(c_int) :: dim = 0
    integer(c_int) :: order = 0
    real(c_double) :: tol = 0
    real(c_double) :: minimum_leaf_fraction = 0
    integer(c_int) :: split_multi_eval = 0
    integer(c_int) :: max_depth = 50
  end type baobzi_input_t

  interface
    function baobzi_init (input, center, half_length) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      import :: baobzi_input_t
      type(baobzi_input_t), intent(in) :: input
      real(kind=c_double), dimension(*) :: center, half_length
      type(c_ptr) :: func
    end function baobzi_init

    subroutine baobzi_save (fin, fname) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: fin
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_ptr) :: func
    end subroutine baobzi_save

    function baobzi_restore (fname) bind(c) result(func)
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_ptr) :: func
    end function baobzi_restore

    function baobzi_eval (f, x) bind(c) result(y)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: f
      real(kind=c_double), dimension(*) :: x
      real(kind=c_double) :: y
    end function baobzi_eval

    subroutine baobzi_eval_multi (f, x, res, ntrg) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: f
      real(kind=c_double), dimension(*) :: x
      real(kind=c_double), dimension(*) :: res
      integer(kind=c_int), value :: ntrg
    end subroutine baobzi_eval_multi

    subroutine baobzi_stats (fin) bind(c)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: fin
      type(c_ptr) :: func
    end subroutine baobzi_stats

    function baobzi_free (f) bind(c) result(f_null)
      use, intrinsic :: iso_c_binding
      type(c_ptr), intent(in), value :: f
      type(c_ptr) :: f_null
    end function baobzi_free
  end interface

end module baobzi
