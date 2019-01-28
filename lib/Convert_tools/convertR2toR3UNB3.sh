#!/bin/bash
for f in $(ls *.rnx)
do
    ./gfzrnx -f -finp $f -fout ::INP:: --version_out 2
done
