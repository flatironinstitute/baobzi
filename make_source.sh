#!/usr/bin/env bash

for iset in {0..2}; do
    for dim in {1..3}; do
        for order in {6..16..2}; do
            srcfile=src/baobzi_${dim}d_${order}_${iset}.cpp
            echo $srcfile
            printf '#include "baobzi.h"\n' > $srcfile
            printf '#include "baobzi/macros.h"\n' >> $srcfile
            printf '#include "baobzi.hpp"\n\n'  >> $srcfile
            printf 'extern "C" {\n' >> $srcfile
            printf 'BAOBZI_DEFS(%d, %d, %d)\n' $dim $order $iset >> $srcfile
            printf '}\n' >> $srcfile
        done
    done
done
