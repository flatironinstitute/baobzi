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

* [Example use cases](#example-use-cases)
* [Limitations](#limitations)
* [Features](#features)
* [Quick install](#quick-install)
* [Building/testing](#building-testing)
* [Input parameters](#input-parameters)
* [Running with...](#running-with)
    + [C](#c)
    + [C++](#c-1)
    + [Python](#python)
    + [Julia](#julia)
    + [MATLAB](#matlab)
    + [Fortran](#fortran)
* [Environment](#environment)
* [Including in your CMake project](#including-in-your-cmake-project)
* [Roadmap](#roadmap)
* [Known Issues. IMPORTANT PLEASE READ](#known-issues-important-please-read)
* [Why the name?](#why-the-name)

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
  optional static library supported for building C/C++ codes where you don't want to load in
  the shared `baobzi` object, but would rather throw it right in your binary. See [Including in
  your CMake project](#including-in-your-cmake-project).

## Quick install
### Python/Conda
```bash
# create a virtualenv if you want
# I advise --system-site-packages when file quotas can a problem
python3 -m venv --system-site-packages myenv
source myenv/bin/activate
pip install baobzi
```

Or with conda
```bash
# create a conda environment, if you want
conda create -y -n baobzi
conda activate baobzi
pip install baobzi
```

## Building/testing/other languages
Baobzi's only dependencies are cmake >= 3.14, and a C/C++17 compiler (gcc only really,
currently). I get the best performance out of g++-11 right now. While there is a header only
library for C++, it can be quite finicky. Therefore, for optimal performance, I currently
suggest using the C shared/static library interface with gcc rather than the C++ header
directly. See [examples/C/baobzi_timing.c](examples/C/baobzi_timing.c) for an example. On my
Xeon Gold 6128 using one core, this example gets roughly 20M evals/s on a simple 2D example,
and 3M evals/s on a simple 3D example.

```bash
# At FI -- module load gcc cmake matlab
export BAOBZI_ROOT=$HOME/local/baobzi
git clone --recursive https://github.com/flatironinstitute/baobzi.git
cd baobzi
mkdir build
cd build
# Don't supply -march!!! Builds with CPU dispatch
cmake -DBAOBZI_BUILD_MATLAB=True -DCMAKE_INSTALL_PREFIX=$BAOBZI_ROOT ..
make -j $((2*$(nproc)))
make install
```

## Input parameters
Baobzi only has a few input parameters, but they can greatly impact the performance and are
worth playing with for your specific function.

* `func`: Function you want to be approximated
* `dim`: Number of independent variables to your function.
* `order`: Polynomial order used to represent a chunk of your function. Higher order is slower,
  especially in higher dimensions. An evaluation, ignoring search/cache issues, takes
  `O(ORDER^DIM)` time. Search isn't free though, and `baobzi` typically needs fewer
  subdivisions for higher orders, so your function might be faster _and_ use less memory if
  you use a higher order.
* `data`: This parameter is only relevant to C/C++/Fortran. If the function you're fitting
  takes parameters, pack that somehow, and `data` is simply a pointer to that packed
  info. See examples.
* `tol`: The maximum desired relative error between your function and the approximant. It is
  impossible to guarantee that all function evaluations will meet this tolerance, so it's
  important to test the results on the domain you're interested in to ensure that results are to your satisfaction.
* `minimum_leaf_fraction`: Baobzi internally is represented by a tree. However, to speed up
  tree lookups, that tree is divided into subtrees that start at some depth. That depth, by
  default, is one above the first level to have a "leaf." A leaf is a terminal box where the
  function evaluation happens (other nodes just contain pointers to their children). This
  scheme can adversely impact performance in some cases though (such as cases where only one
  node on a level is a leaf). This parameter sets a requirement that baobzi keeps subdividing
  entire levels when the fraction of leaves on a given level is less than this threshold.

  An easy way to think about this is if the parameter is `0.0`, then baobzi will never make a
  leaf node a parent node, and the tree will be as small as possible (but not necessarily well
  balanced). If the parameter is `1.0`, then baobzi will ensure that the final depth of the
  tree is entirely filled with leaves. This tree is perfectly balanced, and therefore exactly a
  uniform grid. Anything between will vary between these two extremes. This parameter can
  EXTREMELY impact performance, especially on 1D trees.
* `split_multi_eval`: When evaluating a vector of points, `baobzi` can currently use one of two
  strategies which can dramatically impact performance, depending on the tree and computer. The
  default is to split the evaluation of the points into two stages, one where the boxes are
  calculated in one pass, and then the points are evaluated with them in a second. This tends
  to help when there are a large number of nodes, as it increases the chance of a cache hit by
  not ever loading in extra evaluation data. However it requires making a temporary data
  structure to hold this, which costs time/memory.

  The second is to just brute force evaluate the points as they appear in the target
  order. This typically works very well in 1D, for small trees, but otherwise has poor performance.

  Setting this parameter to 1 uses the typically faster 'split' model, while setting it to 0
  will use the direct model.
* `max_depth`: Maximum depth allowed for tree before interpolation considered failed.

## Running with...
All examples require your project know where the `baobzi` shared object is located. In the
example above, it's located in either the `$BAOBZI_ROOT/lib` or `$BAOBZI_ROOT/lib64` directory,
depending on your system.
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$BAOBZI_ROOT/lib
export LIBRARY_PATH=$LIBRARY_PATH:$BAOBZI_ROOT/lib
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$BAOBZI_ROOT/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$BAOBZI_ROOT/include
export PYTHONPATH=$PYTHONPATH:$BAOBZI_ROOT/share/baobzi/python
export JULIA_LOAD_PATH=JULIA_LOAD_PATH:$BAOBZI_ROOT/share/baobzi/julia
export MATLABPATH=$MATLABPATH:$BAOBZI_ROOT/share/baobzi/matlab
```

### C
A more complicated example than below: [examples/C/baobzi_timing.c](examples/C/baobzi_timing.c)

```C
// test_baobzi.c
#include <baobzi.h>
#include <stdio.h>

double testfunc(const double *x, const void *data) { return x[0] * x[1]; }

int main(int argc, char *argv[]) {
    baobzi_input_t input = {
        .func = testfunc,
        .data = NULL,
        .dim = 2,
        .order = 6,
        .tol = 1E-10,
        .minimum_leaf_fraction = 0.0,
        .split_multi_eval = 0,
        .max_depth = 50
    };

    const double hl[2] = {1.0, 1.0};
    const double center[2] = {0.0, 0.0};
    const double x[2] = {0.25, 0.25};

    baobzi_t func_approx = baobzi_init(&input, center, hl);
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

### C++
A more complicated example than below: [examples/c++/baobzi_timing.cpp](examples/c++/baobzi_timing.cpp)
```c++
// test_baobzi.cpp
#include <baobzi.hpp>
#include <cstdio>

double testfunc(const double *x) { return x[0] * x[1]; }

int main(int argc, char *argv[]) {
    using baobzi::Baobzi;
    baobzi_input_t input = {
        .func = testfunc,
        .data = NULL,
        .dim = 2,
        .order = 6,
        .tol = 1E-10,
        .minimum_leaf_fraction = 0.0,
        .split_multi_eval = 0,
        .max_depth = 50
    };

    const double hl[2] = {1.0, 1.0};
    const double center[2] = {0.0, 0.0};
    const double x[2] = {0.25, 0.25};
    {
        Baobzi func_approx(&input, center, hl);
        printf("%g\n", func_approx(x));
        func_approx.save("func_approx.baobzi");
    }

    Baobzi func_approx("func_approx.baobzi");
    printf("%g\n", func_approx(x));

    return 0;
}
```

```bash
g++ -o test_baobzi.cpp -lbaobzi
```

### Python
[examples/python/simple2d.py](examples/python/simple2d.py)
```python3
# simple2d.py
from baobzi import Baobzi
import numpy as np

def py_test_func(x):
    return x[0] * x[1]

center = np.array([0.0, 0.0])
hl = np.array([1.0, 1.0])
point = np.array([0.25, 0.25])
tol = 1E-8
minimum_leaf_fraction = 0.0 # optional/default
split_multi_eval = 1 # optional/default
max_depth = 50 # optional/default

test = Baobzi(py_test_func, 2, 6, center, hl, 1E-8, minimum_leaf_fraction, split_multi_eval, max_depth)
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
[examples/julia/simple2d.jl](examples/julia/simple2d.jl)
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
minimum_leaf_fraction = 0.0 # optional/default
split_multi_eval = 1 # optional/default
max_depth = 50 # optional/default
output_file = "simple2d.baobzi"

func_approx = baobzi.init(testfunc, dim, order, center, hl, tol, minimum_leaf_fraction, split_multi_eval, max_depth)
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
function (in its own `myfunc.m` file).

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
minimum_leaf_fraction = 0.0;
split_multi_eval = 1;
max_depth = 50;

func_approx = baobzi('new', 'testfun', dim, order, center, hl, tol, minimum_leaf_fraction, split_multi_eval, max_depth);
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
[examples/fortran/fortran_example.f90](examples/fortran/fortran_example.f90)
```fortran
program main
  use baobzi
  implicit none
  real(kind=c_double) :: center(2), half_length(2), tol
  real(kind=c_double) :: x(2)
  real(kind=c_double), target :: scale_factor
  type(c_ptr) :: func_approx
  character(len=64) :: fname
  type(baobzi_input_t) :: input

  input%func = c_funloc(testfun)
  input%data = c_loc(scale_factor)
  input%dim = 2
  input%order = 6
  input%tol = 1E-8
  input%minimum_leaf_fraction = 0.0
  input%split_multi_eval = 0
  input%max_depth = 50

  center(:) = 0.0
  half_length(:) = 1.0
  x(:) = 0.25

  fname = trim(adjustl('fortran.baobzi'))//char(0)

  func_approx = baobzi_init(input, center, half_length)
  print *, baobzi_eval(func_approx, x) - testfun(x, scale_factor)

  call baobzi_save(func_approx, fname)
  func_approx = baobzi_free(func_approx)

  func_approx = baobzi_restore(fname)
  print *, baobzi_eval(func_approx, x) - testfun(x, scale_factor)

contains
  function testfun (x, scale_factor) bind(c) result(y)
    use, intrinsic :: iso_c_binding
    implicit none
    real(kind=c_double), dimension(*) :: x
    real(kind=c_double) :: scale_factor
    real(kind=c_double) :: y

    y = scale_factor * x(1) * x(2)
  end function testfun

end program main
```

```bash
gfortran -o fortran_example -I$BAOBZI_ROOT/share/baobzi/fortran fortran_example.f90 -lbaobzi
```

## Environment
Most control flow of Baobzi is handled through the input structure, but there is a single
environment variable: `BAOBZI_ARCH`.  This environment variable lets you manually set the
instruction set used by the evaluation. There is usually no reason to use this, but in some
cases, AVX2 will outperform AVX512, or you might want to do some testing on the impact of the
instruction set. Valid values are (case-insensitive) `GENERIC`, `AVX`, `AVX2`, and `AVX512`.


## Including in your CMake project
Here I've added a git submodule in extern/baobzi, and I build and link in the static library.
You can also set the shared library on, and the static off, though you'll want to ensure the
shared library gets installed with your project. Probably better to use the static. In this
example, it added 23MB to my test executable size.

```cmake
set(BAOBZI_BUILD_SHARED OFF CACHE BOOL "")
set(BAOBZI_BUILD_STATIC ON CACHE BOOL "")
set(BAOBZI_BUILD_EXAMPLES OFF CACHE BOOL "")
set(BAOBZI_BUILD_TESTS OFF CACHE BOOL "")
add_subdirectory(extern/baobzi)

add_executable(baobzi_test src/test.cpp)
target_include_directories(baobzi_test PUBLIC extern/baobzi/include)
target_link_libraries(baobzi_test PUBLIC baobzi_static)
```

## Roadmap
See the [issues](https://github.com/blackwer/baobzi/issues) or [project tracker](https://github.com/blackwer/baobzi/projects/1).

## Known Issues. IMPORTANT PLEASE READ
* `baobzi` dimensions are defined on the *semi-open* interval `[x0, x1)`. Baobzi functions
  will return NAN outside the this interval.
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
* Doesn't do 4+ dimensions (though possible, I haven't worked out the fitting procedures
  yet). Probably don't want to go above 5 dimensions, since each node takes `O(8 * ORDER^D)`
  bytes of memory.

## Why the name?
It's a cute version of baobab, or the tree of life, which is already a very popular project
name. The baobab lives an extraordinarily long time, which this interpolator is designed to
do. Plant it (it's a tree!), and use it again and again. That's about the extent of the
metaphor -- try not to think about it too much.
