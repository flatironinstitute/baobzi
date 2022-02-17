#!/usr/bin/env bash

mkdir -p include/baobzi

casefile=include/baobzi/baobzi_cases.h
casefile_restore=include/baobzi/baobzi_cases_restore.h
declfile=include/baobzi/baobzi_decls.h

rm -f $casefile $declfile $casefile_restore

for iset in {0..3}; do
    for dim in {1..3}; do
        for order in {6..16..2}; do
            srcfile=src/baobzi_${dim}d_${order}_${iset}.cpp
            echo $srcfile
            printf '#include "baobzi_template.hpp"\n' > $srcfile
            printf '#include "baobzi.h"\n' >> $srcfile
            printf '#include "baobzi/macros.h"\n\n' >> $srcfile
            printf "namespace baobzi {\n" >> $srcfile
            printf "template\n" >> $srcfile
            printf "typename Function<%d, %d, %d>::CoeffVec Function<%d, %d, %d>::cosarray_;\n" $dim $order $iset $dim $order $iset >> $srcfile
            printf "template\n" >> $srcfile
            printf "Eigen::PartialPivLU<typename Function<%d, %d, %d>::VanderMat> Function<%d, %d, %d>::VLU_;\n" $dim $order $iset $dim $order $iset>> $srcfile
            printf "}\n\n" >> $srcfile
            printf 'extern "C" {\n' >> $srcfile
            printf 'BAOBZI_DEFS(%d, %d, %d)\n' $dim $order $iset >> $srcfile
            printf '}\n' >> $srcfile

            printf 'BAOBZI_CASE(%d, %d, %d)\n' $dim $order $iset >> $casefile
            printf 'BAOBZI_CASE_RESTORE(%d, %d, %d)\n' $dim $order $iset >> $casefile_restore
            printf 'BAOBZI_DECLS(%d, %d, %d)\n' $dim $order $iset >> $declfile
        done
    done
done
