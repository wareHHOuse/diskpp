#!/bin/sh

# This script runs the diffusion problem on quad and hex meshes to obtain
# the data for the guys from poland. It expects that the diffusion problem
# is run on the unit quad and unit hex with load function 1 and dirichlet
# boundary conditions 0 everywhere. Moreover, it expects as output the files
# matrix.txt and rhs.txt

DIFFUSION=../../../diskpp-build/apps/diffusion/diffusion
MUMPS_DRIVER=./mumps/solvercost
MESHDIR=../../meshes
DATADIR=data

function make_matrices() {
    for deg in $(seq 0 4); do
        for sz in 2 4 8 16 32; do
            $DIFFUSION -k $deg $MESHDIR/2D_quads/diskpp/testmesh-$sz-$sz.quad
            mkdir -p $DATADIR/2D-k$deg
            matfile=$DATADIR/2D-k$deg/matrix-2D-k$deg-$sz-$sz.txt
            rhsfile=$DATADIR/2D-k$deg/rhs-2D-k$deg-$sz-$sz.txt
            mv matrix.txt $matfile
            mv rhs.txt $rhsfile
        done
    done

    for deg in $(seq 0 4); do
        for sz in 2 4 8 16 32; do
            $DIFFUSION -k $deg $MESHDIR/3D_hexa/diskpp/testmesh-$sz-$sz-$sz.hex
            mkdir -p $DATADIR/3D-k$deg
            matfile=$DATADIR/3D-k$deg/matrix-3D-k$deg-$sz-$sz-$sz.txt
            rhsfile=$DATADIR/3D-k$deg/rhs-3D-k$deg-$sz-$sz-$sz.txt
            mv matrix.txt $matfile
            mv rhs.txt $rhsfile
        done
    done
}

function run_mumps() {
    for deg in $(seq 0 4); do
        for sz in 2 4 8 16 32; do
            matfile=$DATADIR/2D-k$deg/matrix-2D-k$deg-$sz-$sz.txt
            rhsfile=$DATADIR/2D-k$deg/rhs-2D-k$deg-$sz-$sz.txt
            $MUMPS_DRIVER $matfile $rhsfile > $DATADIR/2D-k$deg/mumps-2D-k$deg-$sz-$sz-$sz.txt
        done
    done

    for deg in $(seq 0 2); do
        for sz in 2 4 8 16 32; do
            matfile=$DATADIR/3D-k$deg/matrix-3D-k$deg-$sz-$sz-$sz.txt
            rhsfile=$DATADIR/3D-k$deg/rhs-3D-k$deg-$sz-$sz-$sz.txt
            $MUMPS_DRIVER $matfile $rhsfile > $DATADIR/3D-k$deg/mumps-3D-k$deg-$sz-$sz-$sz.txt
        done
    done
}

function package() {
    find $DATADIR -name "matrix*.txt" -exec bzip2 {} \;
    find $DATADIR -name "rhs*.txt" -exec bzip2 {} \;
    find $DATADIR -name "?D-k?" -type d -exec tar -f {}.tar -cv {} \;
}

if [ ! $# -eq 1 ]; then
    echo "Usage: $0 <what>"
    exit
fi

case "$1" in
    "mkmat") make_matrices ;;
    "mumps") run_mumps ;;
    "pkg") package ;;
    *) echo "Unknown action" ;;
esac


