
for deg in $(seq 0 4); do
    for sz in 2 4 8 16 32; do
        ./diffusion -k $deg ../../../diskpp/meshes/3D_hexa/diskpp/testmesh-$sz-$sz-$sz.hex
        mv matrix.txt ../../../diskpp/misc/solvercost/data/matrix-3D-k$deg-$sz-$sz-$sz.txt
        mv rhs.txt ../../../diskpp/misc/solvercost/data/rhs-3D-k$deg-$sz-$sz-$sz.txt
    done
done
