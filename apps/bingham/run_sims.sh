#!/bin/sh

for dir in bingham*; do
    cd $dir
    ./bingham_vector config.lua
    cd ..
done
