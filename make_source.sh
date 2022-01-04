#!/usr/bin/env bash

for dim in {1..3}; do
    for order in {6..16..2}; do
        srcfile=src/baobzi_${dim}d_$order.cpp
        printf '#include "baobzi.h"\n' > $srcfile
        printf '#include "baobzi/macros.h"\n' >> $srcfile
        printf '#include "baobzi.hpp"\n\n'  >> $srcfile
        printf 'extern "C" {\n' >> $srcfile
        printf 'BAOBZI_DEFS(%d, %d)\n' $dim $order >> $srcfile
        printf '}\n' >> $srcfile
    done
done
