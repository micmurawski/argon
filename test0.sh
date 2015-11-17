#!/bin/bash
#echo $1
for i in $(eval echo {$1..$2..$3})
do
echo $i
echo `sed -i '8s/.*/'$i'/' input`
echo `./program input avs.dat output >> dane.txt`
done