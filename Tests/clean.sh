#!/bin/bash

rm -rf ./COMPARISON_VERSION ./SC ./SC_reference

for name in `ls -d test_*`; do
    cd $name/new
    rm -f config.last energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat log SC SC_reference
    cd ../old
    rm -f config.last energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat log SC SC_reference
    cd ../..
done

cd volumeChange
for name in `ls -d ?h ?l`; do
    cd $name/new
    rm -f config.last energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat log SC SC_reference
    cd ../old
    rm -f config.last energy.dat stat.dat movie wl-new.dat cluster_stat.dat cluster.dat log SC SC_reference
    cd ../..
done
cd ..
