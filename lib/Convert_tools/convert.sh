#!/bin/bash
for f in $(ls *.crx)
do
./CRX2RNX $f
done
