# Baobzi
An adaptive fast function approximator based on tree search. Word salad aside, `baobzi` is a
tool to convert very CPU intensive function calculations into relatively cheap ones (at the
cost of memory). This is similar to functions like `chebeval` in `MATLAB`, but can be
significantly faster since the order of the polynomial fit can be much much lower to meet
similar tolerances. It also isn't constrained for use only in `MATLAB`.

Internally, `baobzi` represents your function by a grid of binary/quad/oct/N trees, where the
leaves represent the function in some small sub-box of the function's domain with chebyshev
polynomials. When you evaluate your function at a point with baobzi, it searches the tree for
the box containing your point and evaluates using this approximant.

## Example use cases
* Build complicated or high computational cost function in higher level language. Build a
  `baobzi` function and use it as a drop in replacement to your function. Reap speed benefits.
* Take prior function, export it, and use it in a higher performance language: production or
  prototype.
* Your favorite `scipy` function not supported in `numba` yet? Use `baobzi` instead of porting
  that function.

## Limitations
* Can use a _lot_ of memory... especially on oscillatory or rapidly changing functions. If your
  function is periodic, fit your function on one period. Transform the function if necessary.
* Not suited around singularities. Just don't do it. Use multiple `baobzi` objects around your
  singularity if necessary.
* Doesn't do any bounds checking. *It is up to the user to sanitize input*.

## Features
* Relative language agnosticism. The library has a simple C interface, which is then bound to
  multiple languages (C/C++, Fortran, Python, Julia, MATLAB). This means you can make a heavy duty
  function in MATLAB, interpolate it there, and then call it from another language of your choice.
* CPU dispatch -- baobzi will detect your CPU and run an optimized codepath based on that --
  no intervention by the user.
* No library dependencies. All code necessary to build `baobzi` is in `baobzi`. There is an
  optional static library are supported for building C/C++ codes where you don't want to load
  in the shared `baobzi` object, but would rather throw it right in your binary (though I don't
  recommend this, the static library is _huge_ and will grow with more dimension support).

## Building/testing
Baobzi's only dependencies are cmake >= 3.5, and a C/C++17 compiler (gcc only really,
currently). I get _vastly_ better performance out of g++-11 (*ONLY FOR C++ HEADER ONLY BUILDS*)
right now than any other compiler, and I'm not sure exactly why. It's very finicky, where the
size of the compilation unit matters. Therefore, for optimal performance, currently, I suggest
using the C shared/static library interface with gcc rather than the C++ header directly. See
[examples/baobzi_timing.c](examples/baobzi_timing.c) for an example. On my Xeon Gold 6128 using
one core, this example gets roughly 50M evals/s on a simple 2D example, and 20M evals/s on a
simple 3D example.

```bash
# At FI -- module load gcc cmake matlab
git clone --recursive https://github.com/blackwer/baobzi.git
cd baobzi
mkdir build
cd build
# Don't supply -march!!! Builds with CPU dispatch
cmake -DBUILD_MATLAB=True -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=$HOME/local
make -j $((2*$(nproc)))
make install
```

## Running with...
All examples require your project know where the `baobzi` shared object is located. In the
example above, it's located in `$HOME/local/lib64`.
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/lib64
export LIBRARY_PATH=$LIBRARY_PATH:$HOME/local/lib64
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$HOME/local/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$HOME/local/include
export PYTHONPATH=$PYTHONPATH:$HOME/local/share/baobzi/python
export JULIA_LOAD_PATH=$PYTHONPATH:$HOME/local/share/baobzi/julia
export MATLABPATH=$PYTHONPATH:$HOME/local/share/baobzi/matlab
```

### C/C++ (no C++ bindings yet. Easy to wrap in class though)
```C
// test_baobzi.c
#include <baobzi.h>
#include <stdio.h>

double testfunc(const double *x) { return x[0] * x[1]; }

int main(int argc, char *argv[]) {
    const int dim = 2;
    const int order = 6;
    const double tol = 1E-10;
    const double hl[2] = {1.0, 1.0};
    const double center[2] = {0.0, 0.0};
    const double x[2] = {0.25, 0.25};

    baobzi_t func_approx = baobzi_init(testfunc, dim, order, center, hl, tol);
    printf("%g\n", baobzi_eval(func_approx, x));
    baobzi_save(func_approx, "func_approx.baobzi");
    func_approx = baobzi_free(func_approx);

    func_approx = baobzi_restore("func_approx.baobzi");
    printf("%g\n", baobzi_eval(func_approx, x));
    return 0;
}
```

```bash
gcc -o test_baobzi.c -lbaobzi
./test_baobzi
```

### Python
[python/examples/simple2d.py](python/examples/simple2d.py)
```python3
# simple2d.py
from baobzi import Baobzi

def py_test_func(x):
    return x[0] * x[1]

center = [0.0, 0.0]
hl = [1.0, 1.0]
point = [0.25, 0.25]

test = Baobzi(py_test_func, 2, 6, center, hl, 1E-8)
test.save('test.baobzi')
print(test(point))
del test

test2 = Baobzi(filename='test.baobzi')
print(test2(point))
del test2
```

```bash
python3 simple2d.py
```

### Julia
[julia/examples/simple2d.jl](julia/examples/simple2d.jl)
```julia
# simple2d.jl
import baobzi

function testfunc(xp::Ptr{Float64})::Cdouble
    x = unsafe_load(xp, 1)
    y = unsafe_load(xp, 2)
    return x * y
end

center = [0.0, 0.0]
hl = [0.5, 1.0]
test_point = [0.25, 0.25]
dim = 2
order = 6
tol = 1E-8
output_file = "simple2d.baobzi"

func_approx = baobzi.init(testfunc, dim, order, center, hl, tol)
println(baobzi.eval(func_approx, test_point) - testfunc(pointer(test_point)))

baobzi.save(func_approx, output_file)
baobzi.free(func_approx)

func_approx = baobzi.restore(output_file)
println(baobzi.eval(func_approx, test_point) - testfunc(pointer(test_point)))
```

```bash
julia simple2d.jl
```

### MATLAB
MATLAB initialization does not work for anonymous functions. You must declare an actual
function (in its on `myfunc.m` file).

```matlab
%% testfun.m
function [y] = testfun(x)
  y = x(1) * x(2);
end
```

```matlab
%% simple2d.m
dim = 2;
order = 6;
center = [0.0, 0.0];
hl = [1.0, 1.0];
tol = 1E-8;
func_approx = baobzi('new', 'testfun', dim, order, center, hl, tol);
display(func_approx.eval([0.25, 0.25]))
func_approx.save('simple2d.baobzi');
clear func_approx
func_approx = baobzi('restore', 'simple2d.baobzi');
display(func_approx.eval([0.25, 0.25]))
```

```bash
matlab -batch simple2d
```

### Fortran
[examples/fortran_example.f90](examples/fortran_example.f90)
```fortran
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

  func_approx = baobzi_restore(fname)
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
```

```bash
gfortran -o fortran_example -I$HOME/local/share/baobzi/fortran fortran_example.f90 -lbaobzi
```

## Roadmap
See the [issues](https://github.com/blackwer/baobzi/issues) or [project tracker](https://github.com/blackwer/baobzi/projects/1).

## Known Issues. IMPORTANT PLEASE READ
* No out of bounds handling at all. If you call out of bounds, your program could crash or
  worse. Only you can prevent segfaults (or even worse, wrong answers). Note that `baobzi`
  dimensions are defined on the *semi-open* interval `[x0, x1)`. Calling on the upper
  boundaries will segfault or give wrong answers.
* Baobzi can't handle singularities or otherwise pathological functions like `sin(1/x)`. It
  will eat your memory and crash. I intend to handle these issues, but I want to work on an API
  that allows for more options to deal with them. That will take time.
    * Workaround: for now, just piecewise a bunch of baobzi functions to work around
      singularities. I hope to add an option to add exclusion domain, or a bunch of domains to
      fit and represent as a single function.
* In one and two dimensions, `baobzi` will determine tolerance matching by looking at the
  values of certain chebyshev coefficients. This tends to underestimate the fit, resulting in
  better precision than you expected (though there are exceptions). In 3D, it uses a sampling
  technique. This effectively doubles the time to fit, but will give you a tolerance closer to
  the expected one. Start low, and move up until you find a fit that works for you. Eventually
  this will be an option in the API.
* MATLAB initialization does not work for anonymous functions. You must declare an actual
  function (in its own `myfunc.m` file).

## Why the name?
It's a cute version of baobab, or the tree of life, which is already a very popular project
name. The baobab lives an extraordinarily long time, which this interpolator is designed to
do. Plant it (it's a tree!), and use it again and again. That's about the extent of the
metaphor -- try not to think about it too much.
