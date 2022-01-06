# Baobzi
An adaptive fast function approximator based on tree search.

## Why Baobzi?
It's a cute version of baobab, or the tree of life, which is already a very popular project
name. The baobab lives an extraordinarily long time, which this interpolator is designed to
do. Plant it (it's a tree!), and use it again and again. That's about the extent of the
metaphor -- try not to think about it too much.

## Building/testing
Baobzi's only dependencies are cmake >= 3.10, eigen >= 3.4, and a C++17 compiler. I get
_vastly_ better performance out of g++-11 (ONLY FOR C++ HEADER ONLY BUILDS) right now than any
other compiler, and I'm not sure exactly why. It's very finicky, where the size of the
compilation unit matters. Therefore, for optimal performance, currently, I suggest using the C
shared/static library interface with gcc rather than the C++ header directly. See
[examples/baobzi_timing.c](examples/baobzi_timing.c) for an example. On my Xeon Gold 6128 using
one core, this example gets roughly 50M evals/s on a simple 2D example, and 20M evals/s on a
simple 3D example.

```bash
# At FI -- module load gcc cmake eigen
git clone --recursive https://github.com/blackwer/baobzi.git
cd baobzi
mkdir build
cd build
# -march=broadwell performance is basically the same, and I suggest it for FI resources
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_FLAGS="-march=native" ..
make -j
./baobzi_timing_c
```

## Roadmap
* ~~CPU dispatch. With the code generation in place this shouldn't actually be too hard, I don't
  think. Famous last words.~~
* Fix header-only performance for C++ (or make a C++ non-template API).
* Multiple function evaluation. i.e. fit multiple functions on the same tree, and evaluate and
  return the results from all of them on a single call.
    * Will require some fairly massive refactoring. Benefit could be large, if the tree
      searches are slow, but I don't expect much honestly. Low priority, but worth a look.
* Exclusion volume (don't interpolate in a given region)
    * Will have to carefully think about the best way to do this. Getting the API right is key.
* Add Eigen to submodules. It's header only so _shrug_. Should be a pinch.
* Add serialization. Truly expensive functions should have to only be fit once and re-used.
* Add logging/reporting functionality -- probably with spdlog
    * Memory/tree statistics
* Auto-optimization for 'spiky' trees
* Maximum depth/infinite recursion detection/memory limits
* Fix intel compilation
* Add other language bindings (or direct implementations)
    * Julia
    * Python
    * Fortran
* Add roadmap issues to... the issues tracker
* Allow for choice of sampling or coefficient to judge fit
* Your suggestion here -- I left my hard copy of the roadmap in the office.
