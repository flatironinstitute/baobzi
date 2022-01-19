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
# At FI -- module load gcc cmake
git clone --recursive https://github.com/blackwer/baobzi.git
cd baobzi
mkdir build
cd build
# Don't supply -march!!! Builds with CPU dispatch so avx, avx2, and avx512 are all supported by default
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -j
./baobzi_timing_c
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

## Why the name?
It's a cute version of baobab, or the tree of life, which is already a very popular project
name. The baobab lives an extraordinarily long time, which this interpolator is designed to
do. Plant it (it's a tree!), and use it again and again. That's about the extent of the
metaphor -- try not to think about it too much.
