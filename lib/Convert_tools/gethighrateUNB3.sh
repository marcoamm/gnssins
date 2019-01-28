#!/bin/bash
for d in {306..335}
do
   for h in 0{1..9} {10..23}
        do
        wget ftp://cddis.gsfc.nasa.gov/gnss/data/highrate/2016/$d/16d/$h/UNB3*
        done
done

