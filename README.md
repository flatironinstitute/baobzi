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
# Don't supply -march!!! Builds with CPU dispatch so avx, avx2, and avx512 are all supported by default
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -j
./baobzi_timing_c
```

## Roadmap
See the [issues](https://github.com/blackwer/baobzi/issues) or [project tracker](https://github.com/blackwer/baobzi/projects/1).
