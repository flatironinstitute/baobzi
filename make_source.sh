#!/bin/bash

mkdir -p include/baobzi

casefile=include/baobzi/baobzi_cases.h
casefile_restore=include/baobzi/baobzi_cases_restore.h
declfile=include/baobzi/baobzi_decls.h

rm -f $casefile $declfile $casefile_restore

dims=$(seq 1 3)
isets=$(seq 0 3)
orders=$(seq 6 2 16)

for iset in ${isets[@]}; do
    srcfile=src/baobzi_${iset}.cpp
    echo $srcfile
    printf '#include "baobzi_template.hpp"\n' > $srcfile
    printf '#include "baobzi.h"\n' >> $srcfile
    printf '#include "baobzi/macros.h"\n\n' >> $srcfile
    printf 'namespace baobzi {\n' >> $srcfile

    if (( $iset > 0 )); then
        printf '#ifdef BAOBZI_CPU_DISPATCH\n' >> $casefile_restore
        printf '#ifdef BAOBZI_CPU_DISPATCH\n' >> $casefile
    fi

    for dim in ${dims[@]}; do
        for order in ${orders[@]}; do
            printf "template\n" >> $srcfile
            printf "typename Function<%d, %d, %d>::VecOrderD Function<%d, %d, %d>::cosarray_;\n" $dim $order $iset $dim $order $iset >> $srcfile
            printf "template\n" >> $srcfile
            printf "Eigen::PartialPivLU<typename Function<%d, %d, %d>::VanderMat> Function<%d, %d, %d>::VLU_;\n\n" $dim $order $iset $dim $order $iset>> $srcfile
        done
    done
    printf "}\n\n" >> $srcfile

    printf 'extern "C" {\n' >> $srcfile
    for dim in ${dims[@]}; do
        for order in ${orders[@]}; do
            printf 'BAOBZI_DEFS(%d, %d, %d)\n' $dim $order $iset >> $srcfile

            printf 'BAOBZI_CASE(%d, %d, %d)\n' $dim $order $iset >> $casefile
            printf 'BAOBZI_CASE_RESTORE(%d, %d, %d)\n' $dim $order $iset >> $casefile_restore
            printf 'BAOBZI_DECLS(%d, %d, %d)\n' $dim $order $iset >> $declfile
        done
    done
    printf '}\n' >> $srcfile

    if (( $iset > 0 )); then
        printf '#endif\n' >> $casefile_restore
        printf '#endif\n' >> $casefile
    fi
done
