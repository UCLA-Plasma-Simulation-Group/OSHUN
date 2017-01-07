#!/bin/sh
echo "Enter directory name:"
read DIR
#DIR="test"
echo "Enter number of nodes:"
read NUM_NODES

mkdir ./bin/$DIR

cp ./bin/oshun1d.e ./input/inputdeck ./bin/$DIR/
cd ./bin/$DIR
mpiexec -np $NUM_NODES ./oshun1d.e
