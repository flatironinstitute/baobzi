# Baobzi
An adaptive fast function approximator based on tree search.

## Why Baobzi?
It's a cute version of baobab, or the tree of life, which is already a very popular project
name. The baobab lives an extraordinarily long time, which this interpolator is designed to
do. Plant it (it's a tree!), and use it again and again. That's about the extent of the
metaphor -- try not to think about it too much.

## Building/testing
Baobzi's only dependencies are cmake >= 3.5, and a C/C++17 compiler (gcc only really,
currently). I get _vastly_ better performance out of g++-11 (ONLY FOR C++ HEADER ONLY BUILDS)
right now than any other compiler, and I'm not sure exactly why. It's very finicky, where the
size of the compilation unit matters. Therefore, for optimal performance, currently, I suggest
using the C shared/static library interface with gcc rather than the C++ header directly. See
[examples/baobzi_timing.c](examples/baobzi_timing.c) for an example. On my Xeon Gold 6128 using
one core, this example gets roughly 50M evals/s on a simple 2D example, and 20M evals/s on a
simple 3D example.

```bash
# At FI -- module load gcc cmake eigen
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

