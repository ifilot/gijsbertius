#!/bin/bash

mkdir -v orbitals

for n in {1..5}; do
    for (( l=0; l<=$n-1; l++ )); do
        for (( m=-$l; m<=$l; m++ )); do
            echo "Building $n / $l / $m"
            ../build/gijsbertius -n $n -l $l -m $m -o orbitals/$n$l$m.ply
            ctmconv orbitals/$n$l$m.ply orbitals/$n$l$m.stl
        done
    done
done
