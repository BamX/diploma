#!/bin/bash
./build.sh
for i in {1..15}
do
   /usr/bin/time -f "$i\t%e" ./run.sh $i > /dev/null
done
