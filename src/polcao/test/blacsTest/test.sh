#!/bin/bash

#for i in {1..10}
#do
#  mpirun -n 2 ./blacsSetupTest 16 16
#done
echo "2"
mpirun -n 2 ./blacsSetupTest 64 8
echo "3"
mpirun -n 3 ./blacsSetupTest 64 8
echo "4"
mpirun -n 4 ./blacsSetupTest 64 8
