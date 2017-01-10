#!/bin/sh

echo "2D test cases"
for deg in $(seq 0 4); do
    for sz in 2 4 8 16 32; do
        echo "Running $deg $sz"
        mfile=data/2D-k$deg/matrix-2D-k$deg-$sz-$sz.txt
        rfile=data/2D-k$deg/rhs-2D-k$deg-$sz-$sz.txt
        ./solvercost $mfile $rfile > mumps-2D-k$deg-$sz-$sz.txt
    done
done

echo "3D test cases"
for deg in $(seq 0 2); do
    for sz in 2 4 8 16 32; do
        echo "Running $deg $sz"
        mfile=data/3D-k$deg/matrix-3D-k$deg-$sz-$sz-$sz.txt
        rfile=data/3D-k$deg/rhs-3D-k$deg-$sz-$sz-$sz.txt
        ./solvercost $mfile $rfile > mumps-3D-k$deg-$sz-$sz-$sz.txt
    done
done
