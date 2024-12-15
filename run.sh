#!/bin/bash

export KMP_AFFINITY=granularity=fine,compact,1,0

NUMACTL="numactl -i all"

for transa in n t;
do
  for transb in n t;
  do
    for sz in {128..2048..128};
    do
      ${NUMACTL} ./sgemm.bin $transa $transb $sz $sz $sz 1. 1. 2>&1 | tee -a sgemm_${transa}_${transb}.log
    done
  done
done
